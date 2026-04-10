#!/usr/bin/env python3
import numpy as np
import torch
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, precision_score, f1_score, 
                             recall_score, roc_auc_score, balanced_accuracy_score)
from xgboost import XGBClassifier
import os
import joblib

# ------------------------
# 1. Data Loading Logic (Same as CNN)
# ------------------------
class BinaryClfData:
    def __init__(self, virus_input_ids_path, non_virus_input_ids_path, extra_neg_paths_dict, neg_pos_ratio=5):
        print("Loading main data...")
        # map_location='cpu' prevents GPU memory usage during loading
        self.virus = torch.load(virus_input_ids_path, map_location='cpu')
        self.non_virus = torch.load(non_virus_input_ids_path, map_location='cpu')
        
        # 1. Load three extra negative samples
        self.extra_data = {}
        for name, path in extra_neg_paths_dict.items():
            if os.path.exists(path):
                print(f"Loading extra negative ({name}) from {path}...")
                self.extra_data[name] = torch.load(path, map_location='cpu')
            else:
                print(f"[Warning] Path not found for {name}: {path}")

        # 2. Sample original negative samples (Uniref90)
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

        # ================= Stratified Split =================
        train_inputs = []
        train_labels = []
        val_inputs = []
        val_labels = []

        # Helper function: split and collect
        def split_and_collect(data_tensor, label_value, name):
            if len(data_tensor) == 0:
                print(f"Warning: Dataset {name} is empty, skipping.")
                return 0, 0
            
            # Convert to NumPy arrays for XGBoost
            # XGBoost input is typically (N_samples, N_features)
            # Here input_ids are shaped (N, 1024), used directly as features
            X = data_tensor.numpy()
            y = np.full((X.shape[0],), label_value)
            
            # Split (80% Train, 20% Val)
            X_tr, X_val, y_tr, y_val = train_test_split(
                X, y, test_size=0.2, random_state=42
            )
            
            train_inputs.append(X_tr)
            train_labels.append(y_tr)
            val_inputs.append(X_val)
            val_labels.append(y_val)
            
            print(f"  -> {name}: Train={len(X_tr)}, Val={len(X_val)}")
            return len(X_tr), len(X_val)

        print("\nSplitting datasets individually to maintain balance...")
        
        # 3.1 Positive samples (Virus) -> Label 1
        split_and_collect(self.virus, 1, "Virus (Positive)")
        
        # 3.2 Base negative samples (Uniref90) -> Label 0
        split_and_collect(non_virus_subset, 0, "Uniref90 (Negative)")
        
        # 3.3 Extra negative samples -> Label 0
        for name, data in self.extra_data.items():
            split_and_collect(data, 0, f"Extra_{name} (Negative)")

        # 4. Merge all data
        self.X_train = np.concatenate(train_inputs)
        self.y_train = np.concatenate(train_labels)
        self.X_val = np.concatenate(val_inputs)
        self.y_val = np.concatenate(val_labels)

        print("-" * 30)
        print(f"Final Train Size: {self.X_train.shape}")
        print(f"Final Val Size:   {self.X_val.shape}")
        print(f"Train Class Dist: Pos={np.sum(self.y_train==1)}, Neg={np.sum(self.y_train==0)}")
        print(f"Val Class Dist:   Pos={np.sum(self.y_val==1)}, Neg={np.sum(self.y_val==0)}")
        print("-" * 30)

# ------------------------
# 2. Main Execution
# ------------------------
def main():
    # 1. Data path configuration
    virus_input_ids_path = 'train.pt'
    non_virus_input_ids_path = "input_ids_nonvirus.pt"
    
    extra_neg_paths = {
        "protease": "other_rna_protease_train_data.pt",
        "capsid":   "other_rna_capsid_train_data.pt",
        "helicase": "other_rna_helicase_train_data.pt"
    }

    # 2. Prepare data
    data = BinaryClfData(
        virus_input_ids_path, 
        non_virus_input_ids_path, 
        extra_neg_paths_dict=extra_neg_paths,
        neg_pos_ratio=5
    )

    # 3. Compute scale_pos_weight (for class imbalance)
    # XGBoost recommendation: sum(negative instances) / sum(positive instances)
    num_neg = np.sum(data.y_train == 0)
    num_pos = np.sum(data.y_train == 1)
    scale_pos_weight = num_neg / num_pos
    print(f"Computed scale_pos_weight: {scale_pos_weight:.4f}")

    # 4. Initialize XGBoost model
    print("Initializing XGBoost...")
    # === Key fix ===
    # eval_metric and early_stopping_rounds must be passed here, not in fit()
    model = XGBClassifier(
        n_estimators=2000,          # Maximum number of trees
        learning_rate=0.05,         # Learning rate
        max_depth=8,                # Tree depth
        subsample=0.8,              # Sample subsampling
        colsample_bytree=0.8,       # Feature subsampling
        objective='binary:logistic',
        n_jobs=16,                  # Number of CPU cores
        tree_method='hist',         # Histogram-based acceleration (use 'gpu_hist' for GPU)
        device='cuda' if torch.cuda.is_available() else 'cpu',  # Try to use GPU
        scale_pos_weight=scale_pos_weight,  # Handle imbalance
        eval_metric='logloss',      # Moved here
        early_stopping_rounds=50,   # Moved here
        random_state=42
    )

    # 5. Train
    print("Start Training with Early Stopping...")
    # Do not pass eval_metric in fit()
    model.fit(
        data.X_train, 
        data.y_train,
        eval_set=[(data.X_val, data.y_val)],
        verbose=100  # Print logs every 100 rounds
    )

    # 6. Evaluate
    print("Evaluating on Validation Set...")
    # Predict probabilities
    y_pred_proba = model.predict_proba(data.X_val)[:, 1]
    # Predict classes (default threshold 0.5)
    y_pred = model.predict(data.X_val)

    # Compute metrics
    acc = accuracy_score(data.y_val, y_pred)
    bacc = balanced_accuracy_score(data.y_val, y_pred)
    prec = precision_score(data.y_val, y_pred, zero_division=0)
    rec = recall_score(data.y_val, y_pred, zero_division=0)
    f1 = f1_score(data.y_val, y_pred, zero_division=0)
    auc = roc_auc_score(data.y_val, y_pred_proba)

    print("\n" + "="*30)
    print("Validation Metrics:")
    print(f"Accuracy:          {acc:.4f}")
    print(f"Balanced Accuracy: {bacc:.4f}")
    print(f"Precision:         {prec:.4f}")
    print(f"Recall:            {rec:.4f}")
    print(f"F1 Score:          {f1:.4f}")
    print(f"ROC AUC:           {auc:.4f}")
    print("="*30)

    # 7. Save model
    output_dir = "XGB_checkpoint"
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, "xgb_model.json")  
    model.save_model(save_path)
    print(f"Model saved to {save_path}")

if __name__ == "__main__":
    main()