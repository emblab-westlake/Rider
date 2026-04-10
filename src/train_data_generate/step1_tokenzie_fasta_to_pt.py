import os
import argparse
import torch
from Bio import SeqIO
from datasets import Dataset
from transformers import AutoTokenizer


def load_fasta_sequences(input_path, max_length=1024):
    sequences = []
    max_len = 0
    with open(input_path, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            if len(seq) > max_length:
                seq = seq[:max_length]
            sequences.append(seq)
            max_len = max(max_len, len(seq))
    print(f"Loaded {len(sequences)} sequences, max length = {max_len}")
    return sequences


def tokenize_sequences(sequences, tokenizer, max_length=1024):
    ds = Dataset.from_dict({"data": sequences})

    def tokenize_function(examples):
        return tokenizer(
            examples["data"],
            truncation=True,
            max_length=max_length,
            padding="max_length",
            return_tensors="pt",
        )

    encoded = ds.map(tokenize_function, batched=True, remove_columns=["data"])
    encoded.set_format("pt")
    return encoded["input_ids"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--esm_path", type=str, required=True)
    parser.add_argument("--input_fasta", type=str, required=True)
    parser.add_argument("--output_pt", type=str, required=True)
    parser.add_argument("--max_length", type=int, default=1024)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output_pt), exist_ok=True)

    tokenizer = AutoTokenizer.from_pretrained(args.esm_path)
    sequences = load_fasta_sequences(args.input_fasta, max_length=args.max_length)
    input_ids = tokenize_sequences(sequences, tokenizer, max_length=args.max_length)
    torch.save(input_ids, args.output_pt)
    print(f"Saved tokenized tensor to {args.output_pt}")


if __name__ == "__main__":
    main()