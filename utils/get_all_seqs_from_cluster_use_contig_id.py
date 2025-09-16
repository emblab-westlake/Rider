#!/usr/bin/env python3

import os
from collections import defaultdict

# File paths
mmseq2_file = "/usr/commondata/public/gaoyang/dataset/lucaprot_data/contig_split_1K/all_protein_prodigalgv_mmseq2_cluster_results.tsv"
contig_file = "/root/gaoyang/westlake_emblab/Rider/plot_draw/global_environment/consensus_contigs_ids_rider_>=2.txt"
output_file = "/root/gaoyang/westlake_emblab/Rider/docs/data/matched_proteins_from_contigs.txt"

# Step 1: Build contig → [member_ids] dictionary from MMseq2 file
contig_to_members = defaultdict(list)

with open(mmseq2_file, "r") as infile:
    for line_num, line in enumerate(infile, 1):
        try:
            center_id, member_id = line.strip().split('\t')
            # Extract contig ID from the center_id (everything before the last underscore)
            # Assuming center_id format is: contig_id + "_" + protein_number
            contig_id = "_".join(center_id.split("_")[:-1])
            contig_to_members[contig_id].append(member_id)
        except ValueError:
            print(f"[Line {line_num}] Skipped malformed line.", flush=True)

        if line_num % 1_000_000 == 0:
            print(f"Indexed {line_num} lines...", flush=True)

print(f"✅ Finished indexing MMseq2 file. Total contigs found: {len(contig_to_members)}", flush=True)

# Step 2: Load contig list to match
with open(contig_file, "r") as f:
    target_contigs = set(line.strip() for line in f if line.strip())

print(f"Loaded {len(target_contigs)} target contig IDs.", flush=True)

# Step 3: Collect all member IDs for matching contigs
matched_members = []

for contig_id in target_contigs:
    if contig_id in contig_to_members:
        matched_members.extend(contig_to_members[contig_id])

print(f"✅ Found {len(matched_members)} matched protein IDs.", flush=True)

# Step 4: Write to output
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out:
    for member_id in matched_members:
        out.write(member_id + "\n")

print(f"✅ Output written to: {output_file}", flush=True)