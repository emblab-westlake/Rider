#!/usr/bin/env python3

import os

# Input paths
mmseq2_file = "/usr/commondata/public/gaoyang/dataset/lucaprot_data/contig_split_1K/all_protein_prodigalgv_mmseq2_cluster_results.tsv"
center_ids_file = "/root/gaoyang/westlake_emblab/Rider/plot_draw/global_environment/Rider_gt200_ids.txt"
center_ids_file = "/root/gaoyang/westlake_emblab/Rider/plot_draw/global_environment/consensus_contigs_ids_rider_>=2.txt"  

# Output path
output_file = "/root/gaoyang/westlake_emblab/Rider/docs/data/matched_member_ids_Rider>=2.txt"

# Load target center IDs into a set for fast lookup
with open(center_ids_file, "r") as f:
    target_centers = set(line.strip() for line in f if line.strip())

print(f"Loaded {len(target_centers)} target center IDs.")

# Initialize list to store matched member IDs
matched_members = []

# Read the large MMseq2 cluster file line by line
with open(mmseq2_file, "r") as infile:
    for line_num, line in enumerate(infile, 1):
        try:
            center_id, member_id = line.strip().split('\t')
            # Check if the center ID is in our target list
            if center_id in target_centers:
                matched_members.append(member_id)
        except ValueError:
            # Skip malformed lines
            print(f"[Line {line_num}] Skipped malformed line.")

        # Progress status every 1 million lines
        if line_num % 1_000_000 == 0:
            print(f"Processed {line_num} lines...")

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Write matched member IDs to output file
with open(output_file, "w") as out:
    for member_id in matched_members:
        out.write(member_id + "\n")

print(f"✅ Done! Matched {len(matched_members)} member IDs written to:\n{output_file}")