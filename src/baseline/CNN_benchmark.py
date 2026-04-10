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
    roc_curve
)
import argparse
import numpy as np

# Path configuration
DEFAULT_TEST_DATA_PATH = ""
DEFAULT_WEIGHTS_PATH = ""

# Model definition (CNN)
class CNNModel(nn.Module):
    def __init__(self, num_tokens=33, embedding_dim=64, max_length=1024, num_classes=2, class_weights=None):
        super(CNNModel, self).__init__()
        torch.manual_seed(42)
        self.num_tokens = num_tokens
        self.max_length = max_length
        self.embedding_dim = embedding_dim

        self.embedding = nn.Embedding(num_tokens, embedding_dim)
        self.conv1 = nn.Conv1d(in_channels=embedding_dim, out_channels=128, kernel_size=3, padding=1)
        self.match_channels = nn.Conv1d(in_channels=embedding_dim, out_channels=128, kernel_size=1)

        conv_output_size = max_length * 128

        self.fc1 = nn.Linear(conv_output_size, 64)
        self.fc2 = nn.Linear(64, num_classes)

        if class_weights is not None:
            self.register_buffer('class_weights', torch.tensor(class_weights, dtype=torch.float32))
        else:
            self.class_weights = None

    def forward(self, input_ids, attention_mask=None, labels=None):
        x = self.embedding(input_ids)
        x = x.permute(0, 2, 1)
        residual = self.match_channels(x)
        x = F.relu(self.conv1(x))
        x = x + residual
        x = x.flatten(start_dim=1)
        features = F.relu(self.fc1(x))
        logits = self.fc2(features)

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

def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark CNN Classifier.")
    parser.add_argument("--test_data", type=str, default=DEFAULT_TEST_DATA_PATH, help="Path to the test data .pt file.")
    parser.add_argument("--weights", type=str, default=DEFAULT_WEIGHTS_PATH, help="Path to the model.safetensors file.")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--save_path", type=str, default="evaluation_results_cnn_3.npz", help="Path to save plotting data.")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"Using device: {args.device}")

    print(f"Loading data from {args.test_data}...")
    data = torch.load(args.test_data)
    
    if isinstance(data, dict):
        input_ids = data['input_ids']
        labels = data['labels']
    else:
        print("Error: Invalid data format. Expected a dict with 'input_ids' and 'labels'.")
        return

    input_ids = input_ids.to(args.device).long()
    labels = labels.to(args.device).long()
    
    print(f"Data loaded. Shape: {input_ids.shape}")

    print("Initializing CNN model...")
    model = CNNModel(num_tokens=33, embedding_dim=64, max_length=1024, num_classes=2).to(args.device)
    
    print(f"Loading weights from {args.weights}...")
    try:
        state_dict = load_file(args.weights)
    except:
        state_dict = torch.load(args.weights, map_location='cpu')
    
    new_state_dict = {k.replace("module.", "").replace("model.", ""): v for k, v in state_dict.items()}
    model.load_state_dict(new_state_dict, strict=False)
    model.eval()

    print("Running inference...")
    batch_size = 256
    all_probs_list = []
    all_features_list = []
    
    with torch.no_grad():
        for i in range(0, len(input_ids), batch_size):
            batch_input = input_ids[i : i + batch_size]
            outputs = model(batch_input)
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

    print("\n" + "="*30)
    print("TEST RESULTS (Threshold 0.5)")
    print("="*30)
    
    acc = accuracy_score(y_true, y_pred_default)
    f1 = f1_score(y_true, y_pred_default, zero_division=0)
    precision = precision_score(y_true, y_pred_default, zero_division=0)
    recall = recall_score(y_true, y_pred_default, zero_division=0)
    
    try:
        auc_score = roc_auc_score(y_true, y_prob_pos)
    except ValueError:
        auc_score = 0.0

    print(f"Accuracy : {acc:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall   : {recall:.4f}")
    print(f"F1 Score : {f1:.4f}")
    print(f"AUC ROC  : {auc_score:.4f}")
    
    print("\nConfusion Matrix:")
    cm = confusion_matrix(y_true, y_pred_default)
    print(cm)
    
    print("\nReport:")
    print(classification_report(y_true, y_pred_default, target_names=['Negative', 'Positive'], digits=4))

    print("\n" + "="*30)
    print("SAVING DATA")
    print("="*30)
    
    fpr, tpr, _ = roc_curve(y_true, y_prob_pos)

    print(f"Saving results to {args.save_path}...")
    np.savez(
        args.save_path,
        y_true=y_true,
        y_prob=y_prob_pos,
        y_pred=y_pred_default,
        embeddings=all_features,
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