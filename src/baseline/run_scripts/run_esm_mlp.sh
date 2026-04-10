#!/bin/bash

# Demo-friendly script for running ESM-MLP evaluation on multiple test files
cd /path/to/baseline/scripts

# Output root directory for results
out_root="path/to/save/results"

# Iterate over embedding files
for i in /path/to/test_set/ESM2_35M_embeddings/*.pt; do
  echo "Processing file: $i"

  base="$(basename "$i")"
  stem="${base%.pt}"

  # Extract the set index from the filename, e.g. set1 -> 1
  idx="$(echo "$stem" | sed -n 's/.*_set\([0-9]\+\)_.*/\1/p')"
  idx="${idx:-X}"

  # Create output directory for each set
  out_dir="${out_root}/set${idx}"
  mkdir -p "$out_dir"

  # Output file name
  out_path="${out_dir}/evaluation_result_100_10000_ESM_mlp_${idx}.npz"

  CUDA_VISIBLE_DEVICES=1 python /path/to/esm_mlp_benchmark.py \
    --embeddings_pt "$i" \
    --weights "/path/to/ESM_MLP_simple_mask_checkpoint/checkpoint-25200/model.safetensors" \
    --device cuda \
    --save_path "$out_path"
done