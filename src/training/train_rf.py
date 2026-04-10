#!/usr/bin/env python3
import numpy as np
import torch
import os
import joblib
import time
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier

# --- Path configuration ---
VIRUS_PATH = 'train.pt'
NON_VIRUS_PATH = "input_ids_nonvirus.pt"

# Additional negative sample paths
EXTRA_NEG_PATHS = {
    "protease": "other_rna_protease_train_data.pt",
    "capsid":   "other_rna_capsid_train_data.pt",
    "helicase": "other_rna_helicase_train_data.pt"
}

SAVE_DIR = "rf_checkpoint"
MODEL_NAME = "best_rf_model.joblib"

def load_and_process_data():
    print("Loading main data...")
    # Use map_location='cpu' to avoid GPU-related loading issues
    virus = torch.load(VIRUS_PATH, map_location='cpu')
    non_virus = torch.load(NON_VIRUS_PATH, map_location='cpu')
    
    # 1. Load extra negative samples
    extra_data = {}
    for name, path in EXTRA_NEG_PATHS.items():
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
    
    # 3. Stratified split by source
    # Split each source separately so train and validation sets keep a balanced source distribution
    
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
    print(f"Train Class Dist: Pos={np.sum(y_train==1)}, Neg={np.sum(y_train==0)}")
    
    return X_train, X_val, y_train, y_val

def main():
    # 1. Prepare data
    X_train, X_val, y_train, y_val = load_and_process_data()
    
    # 2. Define model
    print("\nInitializing Random Forest...")
    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=20,
        class_weight='balanced',  # Automatically handle class imbalance
        n_jobs=-1,                # Use all CPU cores
        random_state=42,
        verbose=1
    )
    
    # 3. Train
    print("Start Training...")
    start_time = time.time()
    model.fit(X_train, y_train)
    print(f"Training finished in {time.time() - start_time:.2f}s")
    
    # 4. Evaluate
    print("Predicting on validation set...")
    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]
    
    print("\n--- Evaluation Results ---")
    print(f"Accuracy : {accuracy_score(y_val, y_pred):.4f}")
    print(f"F1 Score : {f1_score(y_val, y_pred):.4f}")
    try:
        print(f"AUC      : {roc_auc_score(y_val, y_prob):.4f}")
    except Exception as e:
        print(f"AUC Error: {e}")
        
    print("Confusion Matrix:")
    print(confusion_matrix(y_val, y_pred))
    
    # 5. Save model
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    
    save_path = os.path.join(SAVE_DIR, MODEL_NAME)
    print(f"\nSaving model to {save_path}...")
    joblib.dump(model, save_path)
    print("Done.")

if __name__ == "__main__":
    main()