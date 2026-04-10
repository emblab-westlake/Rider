#!/bin/bash

# Run CNN evaluation on multiple test files
cd /path/to/baseline/scripts

# Output root directory for results
out_root="path/to/save/results"

# Iterate over token ID files
for i in /path/to/test_set/tokenID/*.pt; do
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
  out_path="${out_dir}/evaluation_result_100_10000_CNN_${idx}.npz"

  CUDA_VISIBLE_DEVICES=0 python /path/to/CNN_benchmark.py \
    --test_data "$i" \
    --weights "/path/to/Baseline_CNN_Extra3Negs_1:5_dynamicw/final_best_model/model.safetensors" \
    --device cuda \
    --save_path "$out_path"
done