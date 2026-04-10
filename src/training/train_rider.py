#!/usr/bin/env python3
import os
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset
from collections import Counter
# CUDA_VISIBLE_DEVICES=0 python /root/gaoyang/Rider/baseline/train_rider_mask.py

# Hugging Face Transformers
from transformers import (
    Trainer,
    TrainingArguments,
    EvalPrediction,
    DefaultDataCollator,
    EarlyStoppingCallback
)
from transformers.modeling_outputs import ModelOutput

# Metrics
from sklearn.metrics import (
    accuracy_score, precision_score, f1_score,
    recall_score, roc_auc_score, balanced_accuracy_score
)

# Environment Configuration
os.environ["WANDB_MODE"] = "offline"

# ------------------------
# Model Definition
# ------------------------
class SimpleClassifier(nn.Module):
    def __init__(self,
                 input_dim=480,
                 output_dim=2,
                 hidden_dim=1420,
                 nhead=8,
                 num_layers=8,
                 dropout=0.3,
                 class_weights=None):
        """
        Transformer-based classifier for embedding sequences.
        """
        super(SimpleClassifier, self).__init__()
        
        # Initialize Transformer Encoder Layer
        try:
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=input_dim,
                nhead=nhead,
                dim_feedforward=hidden_dim,
                dropout=dropout,
                layer_norm_eps=1e-6,
                batch_first=True
            )
            self.batch_first = True
        except TypeError:
            # Fallback for older PyTorch versions
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=input_dim,
                nhead=nhead,
                dim_feedforward=hidden_dim,
                dropout=dropout,
                layer_norm_eps=1e-6
            )
            self.batch_first = False

        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.linear = nn.Linear(input_dim, output_dim)

        # Register class weights as a buffer
        if class_weights is not None:
            cw = torch.tensor(class_weights, dtype=torch.float32)
            self.register_buffer('class_weights', cw)
        else:
            self.class_weights = None

    def forward(self, input_ids=None, labels=None):
        """
        Forward pass of the model.
        """
        x = input_ids

        # 1. Normalize input dimensions
        # If input is [Batch, Dim], treat it as a sequence of length 1: [Batch, 1, Dim]
        # This is compatible with the generated mask data, which is a 2D tensor [N, 480]
        if x.dim() == 2:
            x = x.unsqueeze(1)
        elif x.dim() != 3:
            raise ValueError(f"Unexpected input dimensions: {x.dim()}. Expected 2 or 3.")

        # 2. Transformer encoder forward pass
        if self.batch_first:
            x = self.transformer_encoder(x)
        else:
            x = x.permute(1, 0, 2)
            x = self.transformer_encoder(x)
            x = x.permute(1, 0, 2)

        # 3. Feature pooling
        cls_features = x[:, 0, :]

        # 4. Classification head
        logits = self.linear(cls_features)

        # 5. Loss computation
        loss = None
        if labels is not None:
            if hasattr(self, 'class_weights') and self.class_weights is not None:
                weights = self.class_weights.to(logits.device)
                loss_fct = nn.CrossEntropyLoss(weight=weights)
            else:
                loss_fct = nn.CrossEntropyLoss()
            
            loss = loss_fct(logits, labels.long())
            
            return ModelOutput({'loss': loss, 'logits': logits, 'features': cls_features})

        return ModelOutput({'logits': logits, 'features': cls_features})

# ------------------------
# Dataset and Collator
# ------------------------
class ClfDataset(Dataset):
    def __init__(self, input_ids, labels):
        self.input_ids = input_ids
        self.labels = labels

    def __getitem__(self, index):
        return {'input_ids': self.input_ids[index], 'labels': self.labels[index]}

    def __len__(self):
        return int(self.input_ids.size(0))

class CustomDataCollator(DefaultDataCollator):
    def __call__(self, features):
        input_ids = torch.stack([f['input_ids'] for f in features])
        labels = torch.tensor([int(f['labels']) for f in features], dtype=torch.long)
        return {'input_ids': input_ids, 'labels': labels}

# ------------------------
# Utility: Compute Class Weights
# ------------------------
def compute_class_weights_from_labels(labels, num_classes=2):
    if isinstance(labels, torch.Tensor):
        label_list = [int(x) for x in labels.tolist()]
    else:
        label_list = [int(x) for x in labels]
        
    cnt = Counter(label_list)
    total = sum(cnt.values())
    weights = []
    
    for i in range(num_classes):
        ci = cnt.get(i, 0)
        if ci == 0:
            print(f"Warning: Class {i} has 0 samples. Assigning high weight.")
            w = 100.0
        else:
            w = total / (num_classes * ci)
        weights.append(w)
    
    return np.array(weights, dtype=np.float32)

# ------------------------
# Metrics Calculation
# ------------------------
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
    
    metrics['accuracy'] = accuracy_score(labels, preds)
    metrics['balanced_accuracy'] = balanced_accuracy_score(labels, preds)
    
    for avg in ['macro', 'weighted']:
        metrics[f'precision_{avg}'] = precision_score(labels, preds, average=avg, zero_division=0)
        metrics[f'recall_{avg}'] = recall_score(labels, preds, average=avg, zero_division=0)
        metrics[f'f1_{avg}'] = f1_score(labels, preds, average=avg, zero_division=0)

    try:
        if probs.shape[1] == 2:
            metrics['roc_auc'] = roc_auc_score(labels, probs[:, 1])
        else:
            metrics['roc_auc'] = roc_auc_score(labels, probs, multi_class='ovr', average='macro')
    except ValueError:
        metrics['roc_auc'] = float('nan')

    unique, counts = np.unique(labels, return_counts=True)
    metrics['class_counts'] = {str(int(k)): int(v) for k, v in zip(unique, counts)}

    if probs.shape[1] == 2:
        pos_probs = probs[:, pos_label]
        thresholds = np.linspace(0.0, 1.0, search_grid)
        best_metric_val = -1.0
        best_threshold = 0.5
        
        for thr in thresholds:
            thr_preds = (pos_probs >= thr).astype(int)
            
            if optimize_for == 'f1':
                val = f1_score(labels, thr_preds, average='binary', zero_division=0)
            elif optimize_for == 'precision':
                val = precision_score(labels, thr_preds, average='binary', zero_division=0)
            elif optimize_for == 'recall':
                val = recall_score(labels, thr_preds, average='binary', zero_division=0)
            
            if val > best_metric_val:
                best_metric_val = val
                best_threshold = thr
        
        final_preds = (pos_probs >= best_threshold).astype(int)
        metrics['best_threshold'] = best_threshold
        metrics[f'best_threshold_optim_{optimize_for}'] = best_metric_val
        metrics['thr_precision'] = precision_score(labels, final_preds, average='binary', zero_division=0)
        metrics['thr_recall'] = recall_score(labels, final_preds, average='binary', zero_division=0)
        metrics['thr_f1'] = f1_score(labels, final_preds, average='binary', zero_division=0)

    return metrics

def compute_metrics(p: EvalPrediction):
    return compute_metrics_with_threshold_search(p, optimize_for='f1', pos_label=1)

# ------------------------
# Main Execution
# ------------------------
if __name__ == "__main__":
    # ================= Begin modifications =================
    # 1. Update data paths to the mask dataset directory generated in the previous step
    base_data_dir = "training_data"
    
    train_data_path = os.path.join(base_data_dir, "train.pt")
    train_label_path = os.path.join(base_data_dir, "train_labels.pt")
    validation_data_path = os.path.join(base_data_dir, "validation.pt")
    validation_label_path = os.path.join(base_data_dir, "validation_labels.pt")
   
    print("Loading datasets...")
    try:
        train_dataset_tensor = torch.load(train_data_path)
        train_label_tensor = torch.load(train_label_path)
        validation_dataset_tensor = torch.load(validation_data_path)
        validation_label_tensor = torch.load(validation_label_path)
        
        print(f"Loaded Train Data Shape: {train_dataset_tensor.shape}")
        print(f"Loaded Val Data Shape: {validation_dataset_tensor.shape}")
    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
        print(f"Please check if files exist in {base_data_dir}")
        exit(1)

    # Compute class weights
    computed_weights = compute_class_weights_from_labels(train_label_tensor, num_classes=2)
    print("Computed class weights:", computed_weights)

    # Initialize model (keep Transformer architecture)
    model = SimpleClassifier(
        input_dim=480, 
        output_dim=2, 
        hidden_dim=960,
        nhead=8, 
        num_layers=2, 
        dropout=0.3,
        class_weights=computed_weights
    )

    # Initialize optimizer
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-4)

    # Create datasets
    TrainDataset = ClfDataset(train_dataset_tensor, train_label_tensor)
    ValidationDataset = ClfDataset(validation_dataset_tensor, validation_label_tensor)

    # 2. Update output directory to avoid overwriting previous experiments
    output_dir = "rider_checkpoint"
    training_args = TrainingArguments(
        output_dir=output_dir,
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
        dataloader_num_workers=4, 
        remove_unused_columns=False 
    )
    # ================= End modifications =================

    # Initialize Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=TrainDataset,
        eval_dataset=ValidationDataset,
        data_collator=CustomDataCollator(),
        compute_metrics=compute_metrics,
        optimizers=(optimizer, None), 
        callbacks=[EarlyStoppingCallback(early_stopping_patience=8)]
    )

    # Start training
    print("Starting training...")
    result = trainer.train()
    print("Training completed.")
    print(result)