"""
Script Name: rider_benchmark.py
Description: Load embeddings and evaluate using both the default threshold (0.5) and the optimized threshold.
             Save comprehensive data for plotting confusion matrices, ROC curves, and embeddings for UMAP/t-SNE.
"""

import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from safetensors.torch import load_file
from sklearn.metrics import (
    accuracy_score, 
    precision_score, 
    recall_score, 
    f1_score, 
    roc_auc_score, 
    average_precision_score, 
    confusion_matrix, 
    classification_report,
    precision_recall_curve,
    roc_curve
)
import argparse
import numpy as np

# Default path configuration
DEFAULT_EMBEDDINGS_PATH = ""
DEFAULT_WEIGHTS_PATH = ""

class SimpleClassifier(nn.Module):
    def __init__(self,
                 input_dim=480, 
                 output_dim=2,
                 hidden_dim=960,
                 nhead=8,
                 num_layers=2,
                 dropout=0.3):
        super(SimpleClassifier, self).__init__()
        
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=input_dim,
            nhead=nhead,
            dim_feedforward=hidden_dim,
            dropout=dropout,
            batch_first=True 
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.linear = nn.Linear(input_dim, output_dim)

    def forward(self, x):
        if x.dim() == 2:
            x = x.unsqueeze(1) 
        x = self.transformer_encoder(x)
        features = x[:, 0, :] 
        logits = self.linear(features)
        return logits, features

def find_optimal_threshold(y_true, y_scores):
    """
    Find the best classification threshold based on F1 score.
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    
    best_idx = np.argmax(f1_scores)
    
    if best_idx < len(thresholds):
        best_threshold = thresholds[best_idx]
    else:
        best_threshold = thresholds[-1]
        
    best_f1 = f1_scores[best_idx]
    return best_threshold, best_f1

def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark Classifier V2.")
    parser.add_argument("--embeddings_pt", type=str, default=DEFAULT_EMBEDDINGS_PATH, help="Path to generated embeddings .pt file.")
    parser.add_argument("--weights", type=str, default=DEFAULT_WEIGHTS_PATH, help="Path to model.safetensors.")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--save_path", type=str, default="evaluation_results_rider_1.npz", help="Path to save plotting data.")
    return parser.parse_args()

def load_classifier_weights(model, weight_path):
    if not os.path.exists(weight_path):
        raise FileNotFoundError(f"Weights not found at {weight_path}")
    
    print(f"Loading weights from {weight_path}...")
    try:
        state_dict = load_file(weight_path)
    except Exception:
        state_dict = torch.load(weight_path, map_location='cpu')

    new_state_dict = {}
    for k, v in state_dict.items():
        if k.startswith("model."):
            new_state_dict[k[6:]] = v
        else:
            new_state_dict[k] = v
            
    model.load_state_dict(new_state_dict, strict=False)
    return model

def main():
    args = parse_args()
    
    # 1. Load data
    print(f"Loading test data from {args.embeddings_pt}...")
    data = torch.load(args.embeddings_pt)
    embeddings = data["embeddings"].to(args.device) 
    labels_tensor = data["labels"].to(args.device)
    labels_numpy = data["labels"].numpy()
    
    print(f"Test set size: {len(labels_numpy)}")

    # 2. Initialize model
    model = SimpleClassifier(input_dim=480, output_dim=2).to(args.device)
    model = load_classifier_weights(model, args.weights)
    model.eval()

    # 3. Inference
    print("Running inference...")
    batch_size = 256
    
    all_logits_list = []
    all_probs_list = []
    all_features_list = []
    
    with torch.no_grad():
        for i in range(0, len(embeddings), batch_size):
            batch_emb = embeddings[i : i + batch_size]
            logits, features = model(batch_emb)
            probs = torch.softmax(logits, dim=1)
            
            all_logits_list.append(logits)
            all_probs_list.append(probs)
            all_features_list.append(features.cpu())

    all_logits = torch.cat(all_logits_list, dim=0)
    all_probs = torch.cat(all_probs_list, dim=0)
    all_features = torch.cat(all_features_list, dim=0).numpy()
    print(f"Extracted embeddings shape: {all_features.shape}")

    y_true = labels_numpy
    y_prob_pos = all_probs[:, 1].cpu().numpy()
    y_pred_default = (y_prob_pos >= 0.5).astype(int)

    # 4. Compute metrics for the default threshold
    print("\n" + "="*50)
    print("1. DEFAULT THRESHOLD (0.5) RESULTS")
    print("="*50)

    loss = F.cross_entropy(all_logits, labels_tensor)
    print(f"Test Loss: {loss.item():.4f}")

    try:
        auc_roc = roc_auc_score(y_true, y_prob_pos)
        aupr = average_precision_score(y_true, y_prob_pos)
    except ValueError:
        auc_roc = 0.0
        aupr = 0.0
    
    print(f"AUC-ROC:   {auc_roc:.4f}")
    print(f"AUC-PR:    {aupr:.4f}")
    print(f"Accuracy:  {accuracy_score(y_true, y_pred_default):.4f}")
    print(f"F1 Score:  {f1_score(y_true, y_pred_default, zero_division=0):.4f}")
    print("\nConfusion Matrix (Default):")
    cm_default = confusion_matrix(y_true, y_pred_default)
    print(cm_default)
    print("\nClassification Report (Default):")
    print(classification_report(y_true, y_pred_default, target_names=['Negative', 'Positive'], digits=4))

    # 5. Find the optimal threshold and re-evaluate
    print("\n" + "="*50)
    print("2. OPTIMIZED THRESHOLD RESULTS")
    print("="*50)

    best_thresh, best_f1_theoretical = find_optimal_threshold(y_true, y_prob_pos)
    
    print(f"Found optimal threshold: {best_thresh:.6f}")
    
    y_pred_optimal = (y_prob_pos >= best_thresh).astype(int)
    
    acc_opt = accuracy_score(y_true, y_pred_optimal)
    prec_opt = precision_score(y_true, y_pred_optimal, zero_division=0)
    rec_opt = recall_score(y_true, y_pred_optimal, zero_division=0)
    f1_opt = f1_score(y_true, y_pred_optimal, zero_division=0)

    print(f"Accuracy:  {acc_opt:.4f}")
    print(f"Precision: {prec_opt:.4f}")
    print(f"Recall:    {rec_opt:.4f}")
    print(f"F1 Score:  {f1_opt:.4f}")

    print("\nConfusion Matrix (Optimized):")
    cm_optimal = confusion_matrix(y_true, y_pred_optimal)
    print(cm_optimal)
    
    print("\nClassification Report (Optimized):")
    print(classification_report(y_true, y_pred_optimal, target_names=['Negative', 'Positive'], digits=4))

    # 6. Save data for plotting
    print("\n" + "="*50)
    print("3. SAVING DATA FOR PLOTTING")
    print("="*50)
    
    fpr, tpr, _ = roc_curve(y_true, y_prob_pos)
    
    print(f"Saving results to {args.save_path}...")
    
    np.savez(
        args.save_path,
        y_true=y_true,
        y_prob=y_prob_pos,
        embeddings=all_features,
        y_pred_default=y_pred_default,
        y_pred_optimal=y_pred_optimal,
        metrics_optimal={
            'accuracy': acc_opt,
            'precision': prec_opt,
            'recall': rec_opt,
            'f1': f1_opt,
            'threshold': best_thresh,
            'auc': auc_roc,
            'confusion_matrix': cm_optimal
        },
        roc_curve={
            'fpr': fpr,
            'tpr': tpr,
            'auc': auc_roc
        }
    )
    print("Done. Data saved successfully.")

if __name__ == "__main__":
    main()