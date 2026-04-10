import argparse
import os
import pandas as pd
import torch
from transformers import AutoTokenizer

def parse_args():
    parser = argparse.ArgumentParser(description="Generate token-id dataset from CSV.")
    parser.add_argument(
        "--csv-path",
        type=str,
        default="",
        help="Input CSV file path. Must contain columns: sequence, label"
    )
    parser.add_argument(
        "--output-pt",
        type=str,
        default="",
        help="Output .pt file path"
    )
    parser.add_argument(
        "--model-path",
        type=str,
        default="",
        help="Local tokenizer/model path for AutoTokenizer"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=1024,
        help="Maximum token length for truncation and padding"
    )
    return parser.parse_args()

def generate_token_id_dataset(csv_path: str, output_pt_path: str, model_path: str, max_length: int):
    print(f"Loading tokenizer from {model_path} ...")
    tokenizer = AutoTokenizer.from_pretrained(model_path)

    print(f"Reading CSV from {csv_path} ...")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # Basic cleaning
    original_len = len(df)
    df = df.dropna(subset=['sequence', 'label'])
    print(f"Data loaded. Rows: {len(df)} (Dropped {original_len - len(df)} NaNs)")

    sequences = df['sequence'].astype(str).tolist()
    labels = df['label'].tolist()

    print("Tokenizing sequences...")
    encodings = tokenizer(
        sequences,
        max_length=max_length,
        padding="max_length",
        truncation=True,
        return_tensors="pt"
    )

    input_ids = encodings['input_ids']  # [N, max_length]
    labels_tensor = torch.tensor(labels, dtype=torch.long)

    print(f"Input IDs shape: {input_ids.shape}")
    print(f"Labels shape: {labels_tensor.shape}")

    data_dict = {
        "input_ids": input_ids,
        "labels": labels_tensor
    }

    # Ensure the output directory exists
    out_dir = os.path.dirname(os.path.abspath(output_pt_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"Saving to {output_pt_path} ...")
    torch.save(data_dict, output_pt_path)
    print("Done! Data is ready for testing.")

if __name__ == "__main__":
    args = parse_args()
    generate_token_id_dataset(
        csv_path=args.csv_path,
        output_pt_path=args.output_pt,
        model_path=args.model_path,
        max_length=args.max_length
    )