import os
import random
import argparse

import pandas as pd
from Bio import SeqIO


def load_fasta_to_list(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    seqs = []
    for record in SeqIO.parse(path, "fasta"):
        seqs.append({"id": record.id, "sequence": str(record.seq)})
    return seqs


def save_dataset(data, output_prefix, set_index, pos_n, neg_n):
    random.shuffle(data)

    csv_name = f"{output_prefix}_set{set_index}_p{pos_n}_n{neg_n}.csv"
    fasta_name = f"{output_prefix}_set{set_index}_p{pos_n}_n{neg_n}.fasta"

    df = pd.DataFrame(data)[["id", "sequence", "label"]]
    df.to_csv(csv_name, index=False)

    with open(fasta_name, "w") as f_out:
        for item in data:
            f_out.write(f">{item['id']}\n{item['sequence']}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sample paired benchmark sets from positive and negative FASTA files."
    )
    parser.add_argument("--pos_fasta", required=True, help="Positive FASTA file")
    parser.add_argument("--neg_fasta", required=True, help="Negative FASTA file")
    parser.add_argument("--pos_per_set", type=int, default=100)
    parser.add_argument("--neg_per_set", type=int, default=4000)
    parser.add_argument("--num_sets", type=int, default=3)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()

    total_pos_needed = args.pos_per_set * args.num_sets
    total_neg_needed = args.neg_per_set * args.num_sets

    pos_data = load_fasta_to_list(args.pos_fasta)
    neg_data = load_fasta_to_list(args.neg_fasta)

    if len(pos_data) < total_pos_needed:
        raise ValueError(f"Not enough positive samples: need {total_pos_needed}, got {len(pos_data)}")
    if len(neg_data) < total_neg_needed:
        raise ValueError(f"Not enough negative samples: need {total_neg_needed}, got {len(neg_data)}")

    if args.seed is not None:
        random.seed(args.seed)

    random.shuffle(pos_data)
    random.shuffle(neg_data)

    for i in range(args.num_sets):
        pos_start = i * args.pos_per_set
        pos_end = (i + 1) * args.pos_per_set
        neg_start = i * args.neg_per_set
        neg_end = (i + 1) * args.neg_per_set

        current_pos = pos_data[pos_start:pos_end]
        current_neg = neg_data[neg_start:neg_end]

        for item in current_pos:
            item["label"] = 1
            item["id"] = f"Set{i+1}_POS_{item['id']}"

        for item in current_neg:
            item["label"] = 0
            item["id"] = f"Set{i+1}_NEG_{item['id']}"

        current_set = current_pos + current_neg
        save_dataset(current_set, args.output_prefix, i + 1, args.pos_per_set, args.neg_per_set)

    print("Done.")
    print(f"Unused positive samples: {len(pos_data) - total_pos_needed}")
    print(f"Unused negative samples: {len(neg_data) - total_neg_needed}")


if __name__ == "__main__":
    main()