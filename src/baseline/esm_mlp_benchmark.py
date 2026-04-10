"""
Script Name: esm_mlp_benchmark.py
Description: Evaluate SimpleLinearModel (MLP) on ESM embeddings.
             Save comprehensive data for later plotting of confusion matrices and ROC curves.
"""

import os
import torch
import torch.nn as nn
import torch.nn.functional as F
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

# Model definition
class SimpleLinearModel(nn.Module):
    def __init__(self,
                 input_dim=480,
                 output_dim=2,
                 class_weights=None):
        super(SimpleLinearModel, self).__init__()
        
        self.mlp = nn.Linear(input_dim, output_dim)

        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids=None, labels=None):
        x = input_ids
        
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
            return {'loss': loss, 'logits': logits, 'features': features}
            
        return {'logits': logits, 'features': features}

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
    parser = argparse.ArgumentParser(description="Benchmark MLP Classifier.")
    parser.add_argument("--embeddings_pt", type=str, default=DEFAULT_EMBEDDINGS_PATH, help="Path to the generated embeddings .pt file.")
    parser.add_argument("--weights", type=str, default=DEFAULT_WEIGHTS_PATH, help="Path to model weights (.pt or .safetensors).")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--save_path", type=str, default="evaluation_results_ESM_mlp_2.npz", help="Path to save plotting data.")
    return parser.parse_args()

def load_classifier_weights(model, weight_path):
    if not os.path.exists(weight_path):
        print(f"Warning: Weights not found at {weight_path}. Initializing with random weights (for testing only).")
        return model
    
    print(f"Loading weights from {weight_path}...")
    try:
        state_dict = torch.load(weight_path, map_location='cpu')
        
        if isinstance(state_dict, dict) and 'model_state_dict' in state_dict:
            state_dict = state_dict['model_state_dict']
        elif isinstance(state_dict, dict) and 'state_dict' in state_dict:
            state_dict = state_dict['state_dict']
            
        new_state_dict = {}
        for k, v in state_dict.items():
            if k.startswith("module."):
                new_state_dict[k[7:]] = v
            else:
                new_state_dict[k] = v
                
        model.load_state_dict(new_state_dict, strict=False)
        print("Weights loaded successfully.")
    except Exception as e:
        print(f"Error loading weights: {e}")
        
    return model

def main():
    args = parse_args()
    
    # Load data
    print(f"Loading test data from {args.embeddings_pt}...")
    try:
        data = torch.load(args.embeddings_pt)
        embeddings = data["embeddings"].to(args.device) 
        labels_tensor = data["labels"].to(args.device)
        labels_numpy = data["labels"].numpy()
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    print(f"Test set size: {len(labels_numpy)}")
    print(f"Embedding shape: {embeddings.shape}")

    # Initialize model
    input_dim = embeddings.shape[-1]
    print(f"Initializing MLP with input_dim={input_dim}...")
    
    model = SimpleLinearModel(input_dim=input_dim, output_dim=2).to(args.device)
    model = load_classifier_weights(model, args.weights)
    model.eval()

    # Inference
    print("Running inference...")
    batch_size = 256
    
    all_logits_list = []
    all_probs_list = []
    all_features_list = []
    
    with torch.no_grad():
        for i in range(0, len(embeddings), batch_size):
            batch_emb = embeddings[i : i + batch_size]
            outputs = model(input_ids=batch_emb) 
            logits = outputs['logits']
            features = outputs['features']
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

    # Default threshold results
    print("\n" + "="*50)
    print("1. DEFAULT THRESHOLD (0.5) RESULTS (MLP)")
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

    # Optimized threshold results
    print("\n" + "="*50)
    print("2. OPTIMIZED THRESHOLD RESULTS (MLP)")
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

    # Save data for plotting
    print("\n" + "="*50)
    print("3. SAVING DATA FOR PLOTTING")
    print("="*50)
    
    fpr, tpr, _ = roc_curve(y_true, y_prob_pos)
    
    print(f"Saving results to {args.save_path}...")
    
    np.savez(
        args.save_path,
        y_true=y_true,
        y_prob=y_prob_pos,
        y_pred_default=y_pred_default,
        y_pred_optimal=y_pred_optimal,
        embeddings=all_features,
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