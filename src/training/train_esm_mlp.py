#!/usr/bin/env python3
import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset
import torch.nn.functional as F
from transformers import Trainer, TrainingArguments, EvalPrediction
from sklearn.metrics import (accuracy_score, precision_score, f1_score,
                             recall_score, roc_auc_score, balanced_accuracy_score)
import torch.optim as optim
from transformers import DefaultDataCollator, EarlyStoppingCallback
from collections import Counter
import os

# CUDA_VISIBLE_DEVICES=1 python /root/gaoyang/Rider/baseline/train_esm_mlp_simple.py
# Set wandb to offline mode
os.environ["WANDB_MODE"] = "offline"

# ------------------------
# Model (Simple Linear Only)
# ------------------------
class SimpleLinearModel(nn.Module):
    def __init__(self,
                 input_dim=480,
                 output_dim=2,
                 class_weights=None):
        super(SimpleLinearModel, self).__init__()
        
        # Fix input/output dimensions
        self.mlp = nn.Linear(input_dim, output_dim)

        # Store class weights as a buffer
        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids=None, labels=None):
        x = input_ids
        
        # Handle input shape
        if x.dim() == 2:
            features = x
        elif x.dim() == 3:
            features = x[:, 0, :]
        else:
            raise ValueError(f"Unexpected input dims: {x.dim()}")

        logits = self.mlp(features)

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
# Utility: class weights
# ------------------------
def compute_class_weights_from_labels(labels, num_classes=2, device='cpu'):
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
    w = np.array(weights, dtype=np.float32)
    return w

# ------------------------
# Metrics
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
    try:
        metrics['accuracy'] = accuracy_score(labels, preds)
    except Exception:
        metrics['accuracy'] = float('nan')

    for name, fn in [('precision_macro', lambda l, q: precision_score(l, q, average='macro', zero_division=0)),
                     ('precision_weighted', lambda l, q: precision_score(l, q, average='weighted', zero_division=0)),
                     ('recall_macro', lambda l, q: recall_score(l, q, average='macro', zero_division=0)),
                     ('recall_weighted', lambda l, q: recall_score(l, q, average='weighted', zero_division=0)),
                     ('f1_macro', lambda l, q: f1_score(l, q, average='macro', zero_division=0)),
                     ('f1_weighted', lambda l, q: f1_score(l, q, average='weighted', zero_division=0))]:
        try:
            metrics[name] = float(fn(labels, preds))
        except Exception:
            metrics[name] = float('nan')

    try:
        metrics['balanced_accuracy'] = float(balanced_accuracy_score(labels, preds))
    except Exception:
        metrics['balanced_accuracy'] = float('nan')

    try:
        n_classes = probs.shape[1]
        if n_classes == 2:
            metrics['roc_auc'] = float(roc_auc_score(labels, probs[:, 1]))
        else:
            metrics['roc_auc'] = float(roc_auc_score(labels, probs, multi_class='ovr', average='macro'))
    except Exception:
        metrics['roc_auc'] = float('nan')

    try:
        if probs.shape[1] == 2:
            pos_probs = probs[:, pos_label]
            thresholds = np.linspace(0.0, 1.0, search_grid)
            best_metric = -1.0
            best_threshold = 0.5
            best_stats = {}
            for thr in thresholds:
                thr_preds = (pos_probs >= thr).astype(int)
                if optimize_for == 'f1':
                    val = f1_score(labels, thr_preds, average='binary', zero_division=0)
                elif optimize_for == 'precision':
                    val = precision_score(labels, thr_preds, average='binary', zero_division=0)
                elif optimize_for == 'recall':
                    val = recall_score(labels, thr_preds, average='binary', zero_division=0)
                else:
                    val = f1_score(labels, thr_preds, average='binary', zero_division=0)
                if val > best_metric:
                    best_metric = val
                    best_threshold = float(thr)
                    best_stats = {
                        'thr_precision': float(precision_score(labels, thr_preds, average='binary', zero_division=0)),
                        'thr_recall': float(recall_score(labels, thr_preds, average='binary', zero_division=0)),
                        'thr_f1': float(f1_score(labels, thr_preds, average='binary', zero_division=0))
                    }
            metrics['best_threshold'] = best_threshold
            metrics['best_threshold_optim_' + optimize_for] = float(best_metric)
            metrics.update(best_stats)
    except Exception:
        metrics['best_threshold'] = float('nan')
        metrics['best_threshold_optim_' + optimize_for] = float('nan')

    return metrics

def compute_metrics(p: EvalPrediction):
    return compute_metrics_with_threshold_search(p, optimize_for='f1', pos_label=1, search_grid=101)

# ------------------------
# Main Execution
# ------------------------
if __name__ == "__main__":
    # Path config
    base_data_dir = "training_data"
    
    train_data_path = os.path.join(base_data_dir, "train.pt")
    train_label_path = os.path.join(base_data_dir, "train_labels.pt")
    validation_data_path = os.path.join(base_data_dir, "validation.pt")
    validation_label_path = os.path.join(base_data_dir, "validation_labels.pt")
   
    # Load data
    print("Loading data...")
    train_dataset_tensor = torch.load(train_data_path)
    train_label_tensor = torch.load(train_label_path)
    validation_dataset_tensor = torch.load(validation_data_path)
    validation_label_tensor = torch.load(validation_label_path)

    # Compute class weights
    computed_weights = compute_class_weights_from_labels(train_label_tensor, num_classes=2)
    print("Computed class weights:", computed_weights)

    # Initialize model
    model = SimpleLinearModel(
        input_dim=480, 
        output_dim=2, 
        class_weights=computed_weights
    )
    
    print("Model Structure (Simple Linear):")
    print(model)

    # Optimizer
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-4)

    # Dataset wrappers
    TrainDataset = ClfDataset(train_dataset_tensor, train_label_tensor)
    ValidationDataset = ClfDataset(validation_dataset_tensor, validation_label_tensor)

    # Training arguments
    training_args = TrainingArguments(
        output_dir="ESM_MLP_checkpoint", 
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
    )

    # Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=TrainDataset,
        eval_dataset=ValidationDataset,
        data_collator=CustomDataCollator(),
        compute_metrics=compute_metrics,
        optimizers=(optimizer, None),
        callbacks=[EarlyStoppingCallback(early_stopping_patience=3)]
    )

    # Start training
    print("Starting Simple Linear Training...")
    result = trainer.train()
    print(result)