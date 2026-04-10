"""
Script Name: 01_generate_embeddings.py
Description: Generate ESM-2 embeddings from a CSV file and save them as a .pt file.
"""

import os
import torch
import pandas as pd
from transformers import AutoTokenizer, EsmModel
from tqdm import tqdm
import argparse

'''
# ================= Configuration =================
python 01_generate_embeddings.py \
  --input_csv /path/to/your/input.csv \
  --esm_model_path /path/to/your/esm_model \
  --output_pt ./bench_set3_embeddings.pt
# ===============================================
'''

# Default path configuration
DEFAULT_CSV_PATH = ""
DEFAULT_ESM_PATH = ""
DEFAULT_OUTPUT_PATH = ""

def parse_args():
    parser = argparse.ArgumentParser(description="Generate embeddings from a CSV file.")
    parser.add_argument("--input_csv", type=str, default=DEFAULT_CSV_PATH, help="Path to the input CSV file.")
    parser.add_argument("--esm_model_path", type=str, default=DEFAULT_ESM_PATH, help="Path to the ESM-2 model directory.")
    parser.add_argument("--output_pt", type=str, default=DEFAULT_OUTPUT_PATH, help="Path to save the output .pt file.")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size.")
    parser.add_argument("--max_length", type=int, default=1024, help="Maximum sequence length.")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Load data
    print(f"Loading data from {args.input_csv}...")
    df = pd.read_csv(args.input_csv)
    
    # Basic cleaning to ensure there are no empty sequences
    df = df.dropna(subset=['sequence', 'label'])
    sequences = df['sequence'].tolist()
    labels = df['label'].tolist()
    ids = df['id'].tolist()
    
    print(f"Total samples: {len(sequences)}")

    # 2. Load model and tokenizer
    print(f"Loading ESM-2 model from {args.esm_model_path}...")
    try:
        tokenizer = AutoTokenizer.from_pretrained(args.esm_model_path)
        model = EsmModel.from_pretrained(args.esm_model_path)
    except Exception as e:
        print(f"Error loading model: {e}")
        print("Please make sure the path points to a valid Hugging Face model directory.")
        return

    model.to(args.device)
    model.eval()

    # 3. Inference loop
    all_embeddings = []
    
    print("Starting inference...")
    with torch.no_grad():
        for i in tqdm(range(0, len(sequences), args.batch_size)):
            batch_seqs = sequences[i : i + args.batch_size]
            
            # Tokenize the batch
            encoded = tokenizer(
                batch_seqs, 
                padding=True, 
                truncation=True, 
                max_length=args.max_length, 
                return_tensors="pt"
            )
            
            input_ids = encoded['input_ids'].to(args.device)
            attention_mask = encoded['attention_mask'].to(args.device)
            
            # Forward pass
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            
            # Extract the [CLS] token embedding (first token, index 0)
            # Shape: [Batch, Hidden_Dim]
            cls_embeddings = outputs.last_hidden_state[:, 0, :]
            
            all_embeddings.append(cls_embeddings.cpu())

    # Concatenate all batches
    final_embeddings = torch.cat(all_embeddings, dim=0)
    final_labels = torch.tensor(labels, dtype=torch.long)
    
    # 4. Save results
    data_to_save = {
        "embeddings": final_embeddings,
        "labels": final_labels,
        "ids": ids
    }
    
    if args.output_pt:
        out_dir = os.path.dirname(os.path.abspath(args.output_pt))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        torch.save(data_to_save, args.output_pt)
        print(f"Done! Saved embeddings to {args.output_pt}")
    else:
        print("Error: --output_pt is empty.")
        return

    print(f"Embeddings shape: {final_embeddings.shape}")

if __name__ == "__main__":
    main()