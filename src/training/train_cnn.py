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
# 1. Model Definition (CNN)
# ------------------------
class CNNModel(nn.Module):
    def __init__(self, num_tokens=33, embedding_dim=64, max_length=1024, num_classes=2, class_weights=None):
        super(CNNModel, self).__init__()
        torch.manual_seed(42)
        self.num_tokens = num_tokens
        self.max_length = max_length
        self.embedding_dim = embedding_dim

        # Embedding layer
        self.embedding = nn.Embedding(num_tokens, embedding_dim)

        # Single convolution layer
        self.conv1 = nn.Conv1d(in_channels=embedding_dim, out_channels=128, kernel_size=3, padding=1)

        # Match input channels for residual connection
        self.match_channels = nn.Conv1d(in_channels=embedding_dim, out_channels=128, kernel_size=1)

        # Output size after convolution
        conv_output_size = max_length * 128

        self.fc1 = nn.Linear(conv_output_size, 64)
        self.fc2 = nn.Linear(64, num_classes)

        # Class weights
        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids, attention_mask=None, labels=None):
        # input_ids: [Batch, Seq_Len]
        
        # Convert token ids to embeddings
        x = self.embedding(input_ids)  # [Batch, Seq_Len, Emb_Dim]

        # Reshape for Conv1d: [Batch, Emb_Dim, Seq_Len]
        x = x.permute(0, 2, 1)
        
        # Residual branch
        residual = self.match_channels(x)
        
        # Convolution + activation
        x = F.relu(self.conv1(x))
        
        # Add residual
        x = x + residual
        
        # Flatten for linear layers
        x = x.flatten(start_dim=1)

        # Fully connected layer
        x = F.relu(self.fc1(x))
        
        # Final logits
        logits = self.fc2(x)

        loss = None
        if labels is not None:
            if self.class_weights is not None:
                weights = self.class_weights.to(logits.device)
                loss_fct = nn.CrossEntropyLoss(weight=weights)
            else:
                loss_fct = nn.CrossEntropyLoss()
            
            loss = loss_fct(logits, labels.long())
            return {'loss': loss, 'logits': logits}

        return {'logits': logits}

# ------------------------
# 2. Dataset Logic (with extra negatives)
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
    def __init__(self, virus_input_ids_path, non_virus_input_ids_path, extra_neg_paths_dict, neg_pos_ratio=5):
        print("Loading main data...")
        # Use map_location='cpu' for safe loading on machines without GPU
        self.virus = torch.load(virus_input_ids_path, map_location='cpu')
        self.non_virus = torch.load(non_virus_input_ids_path, map_location='cpu')
        
        # Load extra negative sets
        self.extra_data = {}
        for name, path in extra_neg_paths_dict.items():
            print(f"Loading extra negative ({name}) from {path}...")
            self.extra_data[name] = torch.load(path, map_location='cpu')

        # Sample base negatives (Uniref90)
        num_pos = self.virus.size(0)
        num_neg_target = num_pos * neg_pos_ratio
        
        if self.non_virus.size(0) < num_neg_target:
            print(f"Warning: Not enough Uniref90 samples. Using all {self.non_virus.size(0)}")
            non_virus_subset = self.non_virus
        else:
            torch.manual_seed(42)
            indices = torch.randperm(self.non_virus.size(0))[:num_neg_target]
            non_virus_subset = self.non_virus[indices]

        print(f"Uniref90 sampled size: {non_virus_subset.size(0)}")

        # Stratified split by source
        train_inputs = []
        train_labels = []
        val_inputs = []
        val_labels = []

        def split_and_collect(data_tensor, label_value, name):
            if len(data_tensor) == 0:
                print(f"Warning: Dataset {name} is empty, skipping.")
                return 0, 0
            
            labels = torch.full((data_tensor.size(0),), label_value)
            
            # 80/20 split
            X_tr, X_val, y_tr, y_val = train_test_split(
                data_tensor.numpy(), labels.numpy(), test_size=0.2, random_state=42
            )
            
            train_inputs.append(X_tr)
            train_labels.append(y_tr)
            val_inputs.append(X_val)
            val_labels.append(y_val)
            
            print(f"  -> {name}: Train={len(X_tr)}, Val={len(X_val)}")
            return len(X_tr), len(X_val)

        print("\nSplitting datasets individually to keep balance...")
        
        # Positive samples
        split_and_collect(self.virus, 1, "Virus (Positive)")
        
        # Base negative samples
        split_and_collect(non_virus_subset, 0, "Uniref90 (Negative)")
        
        # Extra negative samples
        for name, data in self.extra_data.items():
            split_and_collect(data, 0, f"Extra_{name} (Negative)")

        # Merge all splits
        X_train = np.concatenate(train_inputs)
        y_train = np.concatenate(train_labels)
        X_val = np.concatenate(val_inputs)
        y_val = np.concatenate(val_labels)

        # Build datasets
        self.TrainDataset = ClfDataset(torch.tensor(X_train), torch.tensor(y_train))
        self.ValDataset = ClfDataset(torch.tensor(X_val), torch.tensor(y_val))
        
        self.class_num = 2
        
        print("-" * 30)
        print(f"Final Train Size: {len(X_train)}")
        print(f"Final Val Size:   {len(X_val)}")
        print(f"Train Class Dist: Pos={np.sum(y_train==1)}, Neg={np.sum(y_train==0)}")
        print(f"Val Class Dist:   Pos={np.sum(y_val==1)}, Neg={np.sum(y_val==0)}")
        print("-" * 30)

# ------------------------
# 3. Utilities (Weights, Metrics, Collator)
# ------------------------
def compute_class_weights_from_labels(labels, num_classes=2):
    """Compute class weights from training labels."""
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
        
        if probs.shape[1] == 2:
            metrics['roc_auc'] = roc_auc_score(labels, probs[:, 1])
        else:
            metrics['roc_auc'] = roc_auc_score(labels, probs, multi_class='ovr', average='macro')
            
        # Threshold search for binary classification
        if probs.shape[1] == 2:
            pos_probs = probs[:, pos_label]
            thresholds = np.linspace(0.0, 1.0, search_grid)
            best_metric = -1.0
            best_threshold = 0.5
            best_stats = {}
            
            for thr in thresholds:
                thr_preds = (pos_probs >= thr).astype(int)
                val = f1_score(labels, thr_preds, average='binary', zero_division=0)
                if val > best_metric:
                    best_metric = val
                    best_threshold = float(thr)
                    best_stats = {
                        'thr_precision': float(precision_score(labels, thr_preds, average='binary', zero_division=0)),
                        'thr_recall': float(recall_score(labels, thr_preds, average='binary', zero_division=0)),
                        'thr_f1': float(val)
                    }
            metrics['best_threshold'] = best_threshold
            metrics['best_threshold_optim_f1'] = float(best_metric)
            metrics.update(best_stats)
            
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
        attention_mask = torch.ones_like(input_ids)
        return {
            'input_ids': input_ids,
            'attention_mask': attention_mask,
            'labels': labels
        }

# ------------------------
# 4. Main Execution
# ------------------------
if __name__ == "__main__":
    # Data paths
    virus_input_ids_path = 'train.pt'
    non_virus_input_ids_path = "input_ids_nonvirus.pt"
    base_data_dir = "training_data"
    # Extra negative paths
    extra_neg_paths = {
        "protease": "other_rna_protease_train_data.pt",
        "capsid":   "other_rna_capsid_train_data.pt",
        "helicase": "other_rna_helicase_train_data.pt"
    }

    # Prepare data
    dataset_wrapper = BinaryClf(
        virus_input_ids_path, 
        non_virus_input_ids_path, 
        extra_neg_paths_dict=extra_neg_paths,
        neg_pos_ratio=5
    )
    
    # Compute class weights
    train_labels = dataset_wrapper.TrainDataset.labels
    computed_weights = compute_class_weights_from_labels(train_labels, num_classes=2)
    print("Computed class weights:", computed_weights)

    # Initialize model
    model = CNNModel(
        num_tokens=33, 
        embedding_dim=64, 
        max_length=1024, 
        num_classes=2, 
        class_weights=computed_weights
    )

    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-4)

    # Training arguments
    training_args = TrainingArguments(
        output_dir="CNN_checkpoint",
        per_device_train_batch_size=256,
        per_device_eval_batch_size=256,
        num_train_epochs=1000,
        logging_strategy='steps',
        logging_steps=400,
        save_strategy='steps',
        save_steps=400,
        eval_strategy="steps",
        eval_steps=400,
        load_best_model_at_end=True,
        metric_for_best_model="best_threshold_optim_f1", 
        greater_is_better=True,
        save_total_limit=5,
    )

    # Initialize Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=dataset_wrapper.TrainDataset,
        eval_dataset=dataset_wrapper.ValDataset,
        data_collator=CustomDataCollator(),
        compute_metrics=compute_metrics,
        optimizers=(optimizer, None),
        callbacks=[EarlyStoppingCallback(early_stopping_patience=10)]  # Higher patience for a harder dataset
    )

    # Start training
    print("Starting training...")
    result = trainer.train()
    print(result)
    
    trainer.save_model(os.path.join(training_args.output_dir, "final_best_model"))