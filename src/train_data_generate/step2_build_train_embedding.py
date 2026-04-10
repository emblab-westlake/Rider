import os
import argparse
import random

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, EsmModel
from transformers.utils import ModelOutput
from safetensors.torch import load_model
from sklearn.model_selection import train_test_split


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


class LRVMForMaskedLM(nn.Module):
    def __init__(self, model_path, output_dim=480):
        super().__init__()
        self.esm = EsmModel.from_pretrained(model_path)
        self.ln = nn.Linear(output_dim, 33)

    def forward(self, input_ids, attention_mask=None):
        output = self.esm(input_ids=input_ids, attention_mask=attention_mask)
        return output.last_hidden_state


class LRVMForClf(nn.Module):
    def __init__(self, esm_model_path, pretrain_path, clf_class=2, output_dim=480):
        super().__init__()
        self.LRVM = LRVMForMaskedLM(esm_model_path, output_dim)
        print(f"Loading weights from {pretrain_path}...")
        load_model(self.LRVM, pretrain_path, strict=False)
        self.ln = nn.Linear(output_dim, clf_class)

    def forward(self, input_ids, attention_mask=None):
        last_hidden_state = self.LRVM(input_ids=input_ids, attention_mask=attention_mask)
        features = last_hidden_state[:, 0, :]
        logits = self.ln(features)
        return ModelOutput({"logits": logits, "features": features})


class ClfDataset(Dataset):
    def __init__(self, input_ids, labels):
        self.input_ids = input_ids.long()
        self.labels = labels.long()

    def __getitem__(self, index):
        return {"input_ids": self.input_ids[index], "labels": self.labels[index]}

    def __len__(self):
        return self.input_ids.size(0)


class BinaryClf:
    def __init__(self, virus_input_ids_path, non_virus_input_ids_path, extra_neg_paths_dict=None, neg_pos_ratio=5, seed=42):
        print("Loading main data...")
        self.virus = torch.load(virus_input_ids_path)
        self.non_virus = torch.load(non_virus_input_ids_path)
        self.extra_neg_paths_dict = extra_neg_paths_dict or {}

        self.extra_data = {}
        for name, path in self.extra_neg_paths_dict.items():
            print(f"Loading extra negative ({name}) from {path}...")
            self.extra_data[name] = torch.load(path)

        num_pos = self.virus.size(0)
        num_neg_target = num_pos * neg_pos_ratio

        if self.non_virus.size(0) < num_neg_target:
            print(f"Warning: Not enough Uniref90 samples. Using all {self.non_virus.size(0)}")
            non_virus_subset = self.non_virus
        else:
            g = torch.Generator()
            g.manual_seed(seed)
            indices = torch.randperm(self.non_virus.size(0), generator=g)[:num_neg_target]
            non_virus_subset = self.non_virus[indices]

        train_inputs, train_labels = [], []
        val_inputs, val_labels = [], []

        def split_and_collect(data_tensor, label_value, name):
            if data_tensor is None or len(data_tensor) == 0:
                print(f"Warning: {name} is empty, skipping.")
                return

            labels = torch.full((data_tensor.size(0),), label_value, dtype=torch.long)
            X_tr, X_val, y_tr, y_val = train_test_split(
                data_tensor.numpy(),
                labels.numpy(),
                test_size=0.2,
                random_state=seed,
                shuffle=True,
            )

            train_inputs.append(X_tr)
            train_labels.append(y_tr)
            val_inputs.append(X_val)
            val_labels.append(y_val)

            print(f"  -> {name}: Train={len(X_tr)}, Val={len(X_val)}")

        print("\nSplitting datasets individually...")
        split_and_collect(self.virus, 1, "Virus (Positive)")
        split_and_collect(non_virus_subset, 0, "Uniref90 (Negative)")

        for name, data in self.extra_data.items():
            split_and_collect(data, 0, f"Extra_{name} (Negative)")

        X_train = np.concatenate(train_inputs, axis=0)
        y_train = np.concatenate(train_labels, axis=0)
        X_val = np.concatenate(val_inputs, axis=0)
        y_val = np.concatenate(val_labels, axis=0)

        self.TrainDataset = ClfDataset(torch.tensor(X_train), torch.tensor(y_train))
        self.ValDataset = ClfDataset(torch.tensor(X_val), torch.tensor(y_val))

        print("-" * 30)
        print(f"Final Train Size: {len(X_train)}")
        print(f"Final Val Size:   {len(X_val)}")
        print(f"Train Class Dist: Pos={int(np.sum(y_train == 1))}, Neg={int(np.sum(y_train == 0))}")
        print(f"Val Class Dist:   Pos={int(np.sum(y_val == 1))}, Neg={int(np.sum(y_val == 0))}")
        print("-" * 30)


def extract_embeddings(model, dataset, tokenizer, device, batch_size=256, num_workers=4, desc="Dataset"):
    print(f"\nStart extracting embeddings for {desc}...")
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)

    all_embeddings = []
    pad_token_id = tokenizer.pad_token_id

    model.eval()
    with torch.no_grad():
        for i, batch in enumerate(loader):
            input_ids = batch["input_ids"].to(device)
            attention_mask = (input_ids != pad_token_id).long().to(device)

            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            features = outputs.features
            all_embeddings.append(features.cpu().numpy())

            if i % 10 == 0:
                print(f"Processed batch {i}/{len(loader)}", end="\r")

    final_embeddings = np.vstack(all_embeddings)
    print(f"\nFinished {desc}. Shape: {final_embeddings.shape}")
    return torch.tensor(final_embeddings)


def save_tensor(tensor, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    torch.save(tensor, path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str, choices=["normal", "hard_negative"], required=True)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--esm_path", type=str, required=True)
    parser.add_argument("--pretrain_path", type=str, required=True)
    parser.add_argument("--virus_input_ids_path", type=str, required=True)
    parser.add_argument("--non_virus_input_ids_path", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--output_prefix", type=str, required=True)
    parser.add_argument("--neg_pos_ratio", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--num_workers", type=int, default=4)
    parser.add_argument("--extra_neg_protease", type=str, default=None)
    parser.add_argument("--extra_neg_capsid", type=str, default=None)
    parser.add_argument("--extra_neg_helicase", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)
    os.makedirs(args.output_dir, exist_ok=True)

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(args.esm_path)

    extra_neg_paths = {}
    if args.mode == "hard_negative":
        if args.extra_neg_protease:
            extra_neg_paths["protease"] = args.extra_neg_protease
        if args.extra_neg_capsid:
            extra_neg_paths["capsid"] = args.extra_neg_capsid
        if args.extra_neg_helicase:
            extra_neg_paths["helicase"] = args.extra_neg_helicase

    dataset = BinaryClf(
        args.virus_input_ids_path,
        args.non_virus_input_ids_path,
        extra_neg_paths_dict=extra_neg_paths if extra_neg_paths else None,
        neg_pos_ratio=args.neg_pos_ratio,
        seed=args.seed,
    )

    print("Initializing model...")
    model = LRVMForClf(args.esm_path, args.pretrain_path, clf_class=2, output_dim=480).to(device)

    embeddings_train = extract_embeddings(
        model,
        dataset.TrainDataset,
        tokenizer,
        device=device,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        desc="Train Set",
    )

    embeddings_val = extract_embeddings(
        model,
        dataset.ValDataset,
        tokenizer,
        device=device,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        desc="Validation Set",
    )

    save_tensor(embeddings_train, os.path.join(args.output_dir, f"{args.output_prefix}_train_embeddings.pt"))
    save_tensor(embeddings_val, os.path.join(args.output_dir, f"{args.output_prefix}_val_embeddings.pt"))
    save_tensor(dataset.TrainDataset.labels, os.path.join(args.output_dir, f"{args.output_prefix}_train_labels.pt"))
    save_tensor(dataset.ValDataset.labels, os.path.join(args.output_dir, f"{args.output_prefix}_val_labels.pt"))

    print("All Done!")


if __name__ == "__main__":
    main()