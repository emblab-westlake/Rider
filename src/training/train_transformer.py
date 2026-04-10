#!/usr/bin/env python3
import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset
import torch.nn.functional as F
from transformers import Trainer, TrainingArguments, EvalPrediction, DefaultDataCollator, EarlyStoppingCallback
from sklearn.metrics import (accuracy_score, precision_score, f1_score,
                             recall_score, roc_auc_score, balanced_accuracy_score)
from sklearn.model_selection import train_test_split
import torch.optim as optim
from collections import Counter
import os

os.environ["WANDB_MODE"] = "offline"

# ------------------------
# 1. Model Definition (Transformer Baseline)
# ------------------------
class TransformerModel(nn.Module):
    def __init__(self, num_tokens=33, embedding_dim=64, max_length=1024, num_classes=2, num_heads=8, num_transformer_layers=1, class_weights=None):
        super(TransformerModel, self).__init__()
        torch.manual_seed(42)
        self.num_tokens = num_tokens
        self.max_length = max_length
        self.embedding_dim = embedding_dim

        # Embedding layer
        self.embedding = nn.Embedding(num_tokens, embedding_dim)

        # Transformer module
        transformer_layer = nn.TransformerEncoderLayer(d_model=embedding_dim, nhead=num_heads, batch_first=False)
        self.transformer_encoder = nn.TransformerEncoder(transformer_layer, num_layers=num_transformer_layers)

        # Fully connected layers
        conv_output_size = max_length * embedding_dim
        self.fc1 = nn.Linear(conv_output_size, 64)
        self.dropout = nn.Dropout(p=0.5)
        self.fc2 = nn.Linear(64, num_classes)

        # Dynamic class weights
        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids, attention_mask=None, labels=None):
        # input_ids: [Batch, Seq_Len]
        
        # 1. Embedding layer
        x = self.embedding(input_ids)  # [Batch, Seq_Len, Emb_Dim]
        
        # 2. Prepare mask
        if attention_mask is not None:
            # PyTorch Transformer expects True for padding (ignore), False for valid tokens
            src_key_padding_mask = (attention_mask == 0)
        else:
            src_key_padding_mask = None

        # 3. Reshape [Batch, Seq, Dim] -> [Seq, Batch, Dim]
        x = x.permute(1, 0, 2) 

        # 4. Transformer block
        x = self.transformer_encoder(x, src_key_padding_mask=src_key_padding_mask)
        
        # 5. Restore original shape
        x = x.permute(1, 0, 2) 

        # Flatten output before the fully connected layers
        x = x.flatten(start_dim=1)

        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        logits = self.fc2(x)

        loss = None
        if labels is not None:
            if self.class_weights is not None:
                weights = self.class_weights.to(logits.device)
                loss_fct = nn.CrossEntropyLoss(weight=weights)
            else:
                loss_fct = nn.CrossEntropyLoss()
            
            loss = loss_fct(logits.view(-1, self.fc2.out_features), labels.view(-1).long())
            return {'loss': loss, 'logits': logits}

        return {'logits': logits}

# ------------------------
# 2. Dataset Logic (Updated with Extra Negatives & Stratified Split)
# ------------------------
class ClfDataset(Dataset):
    def __init__(self, input_ids, labels):
        self.input_ids = input_ids
        self.labels = labels

    def __getitem__(self, index):
        return self.input_ids[index], self.labels[index]
    
    def __len__(self):
        return self.input_ids.size(0)

class BinaryClf:
    def __init__(self, virus_path, non_virus_path, extra_neg_paths):
        print(f"Loading main virus data from: {virus_path}")
        virus = torch.load(virus_path, map_location='cpu')
        
        print(f"Loading main non-virus data from: {non_virus_path}")
        non_virus = torch.load(non_virus_path, map_location='cpu')
        
        # 1. Load extra negative samples
        extra_data = {}
        for name, path in extra_neg_paths.items():
            print(f"Loading extra negative ({name}) from {path}...")
            extra_data[name] = torch.load(path, map_location='cpu')

        # 2. Sample Uniref90 (1:5 ratio)
        torch.manual_seed(42)
        np.random.seed(42)
        
        num_pos = virus.size(0)
        num_neg_target = num_pos * 5  # 1:5 ratio
        
        if non_virus.size(0) < num_neg_target:
            print(f"Warning: Not enough Uniref90 samples. Using all {non_virus.size(0)}")
            non_virus_subset = non_virus
        else:
            indices = torch.randperm(non_virus.size(0))[:num_neg_target]
            non_virus_subset = non_virus[indices]
            
        print("-" * 30)
        print(f"Virus samples: {virus.size(0)}")
        print(f"Uniref90 sampled: {non_virus_subset.size(0)}")
        
        # 3. Stratified split per source
        # Split each data source separately into 80/20, then merge them
        X_train_list = []
        y_train_list = []
        X_val_list = []
        y_val_list = []

        def split_and_collect(tensor_data, label_val, name):
            if len(tensor_data) == 0:
                return
            
            # Convert to numpy
            X = tensor_data.numpy()
            y = np.full(X.shape[0], label_val)
            
            # 80/20 split
            # Note: we cannot use stratify=y here because y contains only one class.
            # Since each source is split independently, the overall split is effectively stratified.
            X_tr, X_v, y_tr, y_v = train_test_split(
                X, y, test_size=0.2, random_state=42
            )
            
            X_train_list.append(X_tr)
            y_train_list.append(y_tr)
            X_val_list.append(X_v)
            y_val_list.append(y_v)
            
            print(f"  -> {name}: Train={len(X_tr)}, Val={len(X_v)}")

        print("\nSplitting datasets individually...")
        
        # (1) Positive samples: Virus -> Label 1
        split_and_collect(virus, 1, "Virus (Pos)")
        
        # (2) Base negative samples: Uniref90 -> Label 0
        split_and_collect(non_virus_subset, 0, "Uniref90 (Neg)")
        
        # (3) Extra negative samples -> Label 0
        for name, data in extra_data.items():
            split_and_collect(data, 0, f"Extra_{name} (Neg)")

        # 4. Merge all data
        X_train = np.concatenate(X_train_list, axis=0)
        y_train = np.concatenate(y_train_list, axis=0)
        X_val = np.concatenate(X_val_list, axis=0)
        y_val = np.concatenate(y_val_list, axis=0)
        
        print("-" * 30)
        print(f"Final Train Size: {X_train.shape}")
        print(f"Final Val Size:   {X_val.shape}")
        
        # 5. Create Dataset objects
        # Convert back to tensors
        self.TrainDataset = ClfDataset(torch.tensor(X_train).long(), torch.tensor(y_train).float())
        self.ValDataset = ClfDataset(torch.tensor(X_val).long(), torch.tensor(y_val).float())
        
        self.class_num = 2
        print('Dataset preparation done.')

# ------------------------
# 3. Utilities
# ------------------------
def compute_class_weights_from_labels(labels, num_classes=2):
    cnt = Counter([int(x) for x in labels.tolist()]) if isinstance(labels, torch.Tensor) else Counter([int(x) for x in labels])
    total = sum(cnt.values())
    weights = []
    for i in range(num_classes):
        ci = cnt.get(i, 0)
        if ci == 0:
            w = float(total)
        else:
            w = total / (num_classes * ci)
        weights.append(w)
    return np.array(weights, dtype=np.float32)

def compute_metrics_with_threshold_search(p: EvalPrediction, optimize_for='f1', pos_label=1, search_grid=100):
    logits = p.predictions
    if isinstance(logits, tuple) or isinstance(logits, list):
        logits = logits[0]
    
    if isinstance(logits, torch.Tensor):
        probs = F.softmax(logits, dim=1).cpu().numpy()
    else:
        probs = F.softmax(torch.tensor(logits), dim=1).numpy()

    labels = p.label_ids
    preds = np.argmax(probs, axis=1)

    metrics = {}
    try:
        metrics['accuracy'] = accuracy_score(labels, preds)
        metrics['balanced_accuracy'] = balanced_accuracy_score(labels, preds)
        metrics['f1_macro'] = f1_score(labels, preds, average='macro', zero_division=0)
        try:
            metrics['roc_auc'] = roc_auc_score(labels, probs[:, 1])
        except:
            metrics['roc_auc'] = 0.0
            
        # Threshold search
        pos_probs = probs[:, pos_label]
        thresholds = np.linspace(0.0, 1.0, search_grid)
        best_metric = -1.0
        best_threshold = 0.5
        
        for thr in thresholds:
            thr_preds = (pos_probs >= thr).astype(int)
            val = f1_score(labels, thr_preds, average='binary', zero_division=0)
            if val > best_metric:
                best_metric = val
                best_threshold = float(thr)
        
        metrics['best_threshold'] = best_threshold
        metrics['best_threshold_optim_f1'] = float(best_metric)
            
    except Exception as e:
        print(f"Metric error: {e}")
        metrics['best_threshold_optim_f1'] = 0.0

    return metrics

def compute_metrics(p: EvalPrediction):
    return compute_metrics_with_threshold_search(p, optimize_for='f1', pos_label=1, search_grid=101)

class CustomDataCollator(DefaultDataCollator):
    def __call__(self, features):
        input_ids = torch.stack([f[0] for f in features]).long()
        labels = torch.tensor([f[1] for f in features]).long()
        
        # 1 means valid token, 0 means padding
        attention_mask = (input_ids != 0).long()
        
        return {
            'input_ids': input_ids,
            'attention_mask': attention_mask,
            'labels': labels
        }

# ------------------------
# 4. Main Execution
# ------------------------
if __name__ == "__main__":
    # 1. Data paths
    virus_input_ids_path = 'train.pt'
    non_virus_input_ids_path = "input_ids_nonvirus.pt"
    
    # Additional negative sample paths
    extra_neg_paths = {
        "protease": "other_rna_protease_train_data.pt",
        "capsid":   "other_rna_capsid_train_data.pt",
        "helicase": "other_rna_helicase_train_data.pt"
    }

    # 2. Prepare data (with extra paths)
    dataset_wrapper = BinaryClf(virus_input_ids_path, non_virus_input_ids_path, extra_neg_paths)
    
    # 3. Compute class weights dynamically
    train_labels = dataset_wrapper.TrainDataset.labels
    computed_weights = compute_class_weights_from_labels(train_labels, num_classes=2)
    print("Computed class weights:", computed_weights)

    # 4. Initialize Transformer model
    model = TransformerModel(
        num_tokens=33, 
        embedding_dim=64, 
        max_length=1024, 
        num_classes=2, 
        num_heads=8, 
        num_transformer_layers=1,
        class_weights=computed_weights
    )

    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-4)

    # 5. Trainer arguments
    # Update output_dir to reflect the new dataset strategy
    output_dir_path = "transformer_checkpoint"
    
    training_args = TrainingArguments(
        output_dir=output_dir_path,
        per_device_train_batch_size=256,
        per_device_eval_batch_size=256,
        num_train_epochs=1000,
        logging_strategy='steps',
        logging_steps=400,
        save_strategy='steps',
        save_steps=400,
        evaluation_strategy="steps",
        eval_steps=400,
        load_best_model_at_end=True,
        metric_for_best_model="best_threshold_optim_f1", 
        greater_is_better=True,
        save_total_limit=3,
        dataloader_num_workers=4  # Increase workers to speed up data loading
    )

    # 6. Initialize Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=dataset_wrapper.TrainDataset,
        eval_dataset=dataset_wrapper.ValDataset,
        data_collator=CustomDataCollator(),
        compute_metrics=compute_metrics,
        optimizers=(optimizer, None),
        callbacks=[EarlyStoppingCallback(early_stopping_patience=5)]  # Slightly increase patience
    )

    # 7. Start training
    print("Starting training (Transformer Baseline with Extra Negatives)...")
    result = trainer.train()
    print(result)
    
    trainer.save_model(os.path.join(training_args.output_dir, "final_best_model"))