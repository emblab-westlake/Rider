#!/usr/bin/env python3
"""
rider_merge_taxonomy.py

Merge foldseek m8 results with a taxonomy table, output the top hit (by tm_score)
for each query and attach the taxonomy information for that top hit.

Usage:
    python rider_merge_taxonomy.py --m8_file /path/to/file.m8 \
        --taxo_file /path/to/taxo.tsv \
        --output_file /path/to/output.tsv

The script expects the taxonomy file to contain an 'Accession' column that can be
matched to the target accession (with optional version suffix removed).
"""

import argparse
import os
import sys
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Merge foldseek m8 top hits with a taxonomy table")
    parser.add_argument("--m8_file", required=False, default="/usr/commondata/public/gaoyang/foldseek_search_dir/rider_gt200_all_mapping_aln.m8",
                        help="Path to the input m8 file")
    parser.add_argument("--taxo_file", required=False, default="/home/gaoyang/Rider/utils/taxo_Rider_train_dat_parsed_taxonomy.tsv",
                        help="Path to the taxonomy file (must contain 'Accession' column)")
    parser.add_argument("--output_file", required=False, default="/usr/commondata/public/gaoyang/foldseek_search_dir/foldseek_top_hit_taxonomy_lucaprot_gt200_all_mapping_aln.tsv",
                        help="Path to the output file")
    return parser.parse_args()

def main():
    args = parse_args()

    m8_file = args.m8_file
    taxo_file = args.taxo_file
    output_file = args.output_file

    # Check that input files exist
    for p in [m8_file, taxo_file]:
        if not os.path.isfile(p):
            print(f"Error: input file not found: {p}", file=sys.stderr)
            sys.exit(2)

    # Step 1: Read m8 file
    m8_cols = [
        "query", "target", "identity", "aln_len", "matches", "mismatches",
        "q_start", "q_end", "t_start", "t_end", "evalue", "tm_score"
    ]
    # Read without header; treat columns as strings initially
    m8_df = pd.read_csv(m8_file, sep="\t", names=m8_cols, dtype=str, header=None)
    # Convert tm_score to numeric for proper sorting (invalid -> NaN)
    m8_df["tm_score"] = pd.to_numeric(m8_df["tm_score"], errors="coerce")

    # For each query, select the hit with the highest tm_score
    top_hits = m8_df.sort_values("tm_score", ascending=False).drop_duplicates("query").copy()

    # Remove version suffix from target (e.g. .1)
    top_hits["target_clean"] = top_hits["target"].str.replace(r"\.\d+$", "", regex=True)

    # Step 2: Read the taxonomy table and strip column names
    taxo = pd.read_csv(taxo_file, sep="\t", dtype=str)
    taxo.columns = taxo.columns.str.strip()

    # Ensure 'Accession' column exists
    if "Accession" not in taxo.columns:
        print("Error: 'Accession' column not found in taxonomy file. Please check the taxonomy file.", file=sys.stderr)
        sys.exit(3)

    # Drop duplicate accessions, keep first occurrence
    taxo = taxo.drop_duplicates(subset="Accession")

    # Step 3: Merge top_hits with taxonomy information by matching cleaned target to Accession
    merged = pd.merge(
        top_hits,
        taxo,
        how="left",
        left_on="target_clean",
        right_on="Accession"
    )

    # Step 4: Prepare output columns: query, target, plus all taxonomy columns except Accession
    other_cols = [c for c in taxo.columns if c != "Accession"]
    output_cols = ["query", "target"] + other_cols
    final_df = merged.reindex(columns=output_cols)

    # Write output (create output directory if necessary)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Output written: {output_file}")

if __name__ == "__main__":
    main()