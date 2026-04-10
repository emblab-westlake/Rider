import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from safetensors.torch import load_file
from sklearn.metrics import (
    accuracy_score, 
    f1_score, 
    classification_report, 
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    average_precision_score
)
import argparse
import numpy as np

# Path configuration
DEFAULT_TEST_DATA_PATH = ""
DEFAULT_WEIGHTS_PATH = ""

# Model definition (Transformer baseline)
class TransformerModel(nn.Module):
    def __init__(self, num_tokens=33, embedding_dim=64, max_length=1024, num_classes=2, num_heads=8, num_transformer_layers=1, class_weights=None):
        super(TransformerModel, self).__init__()
        torch.manual_seed(42)
        self.num_tokens = num_tokens
        self.max_length = max_length
        self.embedding_dim = embedding_dim

        self.embedding = nn.Embedding(num_tokens, embedding_dim)

        transformer_layer = nn.TransformerEncoderLayer(d_model=embedding_dim, nhead=num_heads, batch_first=False)
        self.transformer_encoder = nn.TransformerEncoder(transformer_layer, num_layers=num_transformer_layers)

        conv_output_size = max_length * embedding_dim
        self.fc1 = nn.Linear(conv_output_size, 64)
        self.dropout = nn.Dropout(p=0.5)
        self.fc2 = nn.Linear(64, num_classes)

        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids, attention_mask=None, labels=None):
        x = self.embedding(input_ids)
        
        if attention_mask is not None:
            src_key_padding_mask = (attention_mask == 0)
        else:
            src_key_padding_mask = None

        x = x.permute(1, 0, 2)
        x = self.transformer_encoder(x, src_key_padding_mask=src_key_padding_mask)
        x = x.permute(1, 0, 2)

        x = x.flatten(start_dim=1)
        features = F.relu(self.fc1(x))
        
        x = self.dropout(features)
        logits = self.fc2(x)

        loss = None
        if labels is not None:
            if self.class_weights is not None:
                weights = self.class_weights.to(logits.device)
                loss_fct = nn.CrossEntropyLoss(weight=weights)
            else:
                loss_fct = nn.CrossEntropyLoss()
            
            loss = loss_fct(logits.view(-1, self.fc2.out_features), labels.view(-1).long())
            return {'loss': loss, 'logits': logits, 'features': features}

        return {'logits': logits, 'features': features}

# Helper function: find optimal threshold
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
    parser = argparse.ArgumentParser(description="Benchmark Transformer Classifier.")
    parser.add_argument("--test_data", type=str, default=DEFAULT_TEST_DATA_PATH, help="Path to test data .pt file.")
    parser.add_argument("--weights", type=str, default=DEFAULT_WEIGHTS_PATH, help="Path to model.safetensors.")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--save_path", type=str, default="evaluation_results_transformer_optimal.npz", help="Path to save plotting data.")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"Using device: {args.device}")

    # 1. Load data
    print(f"Loading data from {args.test_data}...")
    data = torch.load(args.test_data)
    
    if isinstance(data, dict):
        input_ids = data['input_ids']
        labels = data['labels']
        
        if 'attention_mask' in data:
            print("Found 'attention_mask' in dataset.")
            attention_mask = data['attention_mask']
        else:
            print("Warning: 'attention_mask' not found. Generating mask from input_ids (assuming 0 is padding).")
            attention_mask = (input_ids != 0).long()
    else:
        print("Error: Data format incorrect.")
        return

    input_ids = input_ids.to(args.device).long()
    labels = labels.to(args.device).long()
    attention_mask = attention_mask.to(args.device).long()
    
    print(f"Data loaded. Input shape: {input_ids.shape}")

    # 2. Initialize model
    print("Initializing Transformer model...")
    model = TransformerModel(
        num_tokens=33, 
        embedding_dim=64, 
        max_length=1024, 
        num_classes=2,
        num_heads=8,
        num_transformer_layers=1
    ).to(args.device)
    
    print(f"Loading weights from {args.weights}...")
    try:
        state_dict = load_file(args.weights)
    except:
        state_dict = torch.load(args.weights, map_location='cpu')
    
    new_state_dict = {k.replace("module.", "").replace("model.", ""): v for k, v in state_dict.items()}
    model.load_state_dict(new_state_dict, strict=False)
    model.eval()

    # 3. Inference
    print("Running inference...")
    batch_size = 256
    all_probs_list = []
    all_features_list = []
    
    with torch.no_grad():
        for i in range(0, len(input_ids), batch_size):
            batch_input = input_ids[i : i + batch_size]
            batch_mask = attention_mask[i : i + batch_size]
            
            outputs = model(batch_input, attention_mask=batch_mask)
            
            logits = outputs['logits']
            features = outputs['features']
            
            probs = torch.softmax(logits, dim=1)
            all_probs_list.append(probs.cpu())
            all_features_list.append(features.cpu())

    all_probs = torch.cat(all_probs_list, dim=0).numpy()
    all_features = torch.cat(all_features_list, dim=0).numpy()
    print(f"Extracted embeddings shape: {all_features.shape}")

    y_true = labels.cpu().numpy()
    y_prob_pos = all_probs[:, 1]
    
    y_pred_default = (y_prob_pos >= 0.5).astype(int)

    # 4. Calculate metrics for the default threshold
    print("\n" + "="*50)
    print("1. DEFAULT THRESHOLD (0.5) RESULTS")
    print("="*50)
    
    acc = accuracy_score(y_true, y_pred_default)
    f1 = f1_score(y_true, y_pred_default, zero_division=0)
    precision = precision_score(y_true, y_pred_default, zero_division=0)
    recall = recall_score(y_true, y_pred_default, zero_division=0)
    
    try:
        auc_score = roc_auc_score(y_true, y_prob_pos)
        aupr = average_precision_score(y_true, y_prob_pos)
    except ValueError:
        auc_score = 0.0
        aupr = 0.0

    print(f"Accuracy : {acc:.4f}")
    print(f"F1 Score : {f1:.4f}")
    print(f"AUC ROC  : {auc_score:.4f}")
    print(f"AUC PR   : {aupr:.4f}")
    
    print("\nConfusion Matrix (Default):")
    cm_default = confusion_matrix(y_true, y_pred_default)
    print(cm_default)
    
    print("\nReport (Default):")
    print(classification_report(y_true, y_pred_default, target_names=['Negative', 'Positive'], digits=4))

    # 5. Find optimal threshold and re-evaluate
    print("\n" + "="*50)
    print("2. OPTIMIZED THRESHOLD RESULTS")
    print("="*50)

    best_thresh, best_f1_theoretical = find_optimal_threshold(y_true, y_prob_pos)
    print(f"Found Optimal Threshold: {best_thresh:.6f}")

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
    
    print("\nReport (Optimized):")
    print(classification_report(y_true, y_pred_optimal, target_names=['Negative', 'Positive'], digits=4))

    # 6. Save data for plotting
    print("\n" + "="*50)
    print("3. SAVING DATA")
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
            'auc': auc_score,
            'confusion_matrix': cm_optimal
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