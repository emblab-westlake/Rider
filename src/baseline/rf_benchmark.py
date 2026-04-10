import os
import torch
import joblib
import argparse
import numpy as np
from sklearn.metrics import (
    accuracy_score, 
    f1_score, 
    classification_report, 
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve
)

# Path configuration
# Use the same test set as the Transformer evaluation to ensure fair comparison
DEFAULT_TEST_DATA_PATH = ""
DEFAULT_MODEL_PATH = ""

def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark Random Forest Classifier.")
    parser.add_argument("--test_data", type=str, default=DEFAULT_TEST_DATA_PATH, help="Path to test data .pt file.")
    parser.add_argument("--model_path", type=str, default=DEFAULT_MODEL_PATH, help="Path to the .joblib model file.")
    parser.add_argument("--save_path", type=str, default="evaluation_results_rf_3.npz", help="Path to save plotting data.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Load data
    print(f"Loading data from {args.test_data}...")
    try:
        data = torch.load(args.test_data)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    if isinstance(data, dict):
        input_ids = data['input_ids']
        labels = data['labels']
    else:
        print("Error: Invalid data format. Expected a dict with 'input_ids' and 'labels'.")
        return

    # Random Forest requires NumPy arrays on CPU
    print("Converting data to NumPy format...")
    X_test = input_ids.cpu().numpy()
    y_true = labels.cpu().numpy()
    
    print(f"Data loaded. Shape: {X_test.shape}")

    # 2. Load model
    print(f"Loading Random Forest model from {args.model_path}...")
    if not os.path.exists(args.model_path):
        print(f"Error: Model file not found at {args.model_path}")
        return
        
    model = joblib.load(args.model_path)
    print("Model loaded successfully.")

    # 3. Inference
    print("Running inference (this may take a moment on CPU)...")
    
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    # 4. Compute metrics
    print("\n" + "="*30)
    print("TEST RESULTS (Random Forest)")
    print("="*30)
    
    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    
    try:
        auc_score = roc_auc_score(y_true, y_prob)
    except ValueError:
        auc_score = 0.0

    print(f"Accuracy : {acc:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall   : {recall:.4f}")
    print(f"F1 Score : {f1:.4f}")
    print(f"AUC ROC  : {auc_score:.4f}")
    
    print("\nConfusion Matrix:")
    cm = confusion_matrix(y_true, y_pred)
    print(cm)
    
    print("\nReport:")
    print(classification_report(y_true, y_pred, target_names=['Negative', 'Positive'], digits=4))

    # 5. Save data for plotting
    print("\n" + "="*30)
    print("SAVING DATA")
    print("="*30)
    
    fpr, tpr, _ = roc_curve(y_true, y_prob)

    print(f"Saving results to {args.save_path}...")
    
    np.savez(
        args.save_path,
        y_true=y_true,
        y_prob=y_prob,
        y_pred=y_pred,
        embeddings=X_test,
        metrics={
            'accuracy': acc,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'auc': auc_score,
            'confusion_matrix': cm
        },
        roc_curve={
            'fpr': fpr,
            'tpr': tpr,
            'auc': auc_score
        }
    )
    print("Done.")

if __name__ == "__main__":
    main()