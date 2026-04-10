#!/bin/bash

# Run Transformer evaluation on multiple test files
cd /path/to/baseline/scripts

# Output root directory for results
out_root="path/to/save/results"

# Iterate over token ID files
for i in /path/to/test_set/tokenID/*.pt; do
  echo "Processing file: $i"

  base="$(basename "$i")"
  stem="${base%.pt}"

  idx="$(echo "$stem" | sed -n 's/.*_set\([0-9]\+\)_.*/\1/p')"
  idx="${idx:-X}"

  out_dir="${out_root}/set${idx}"
  mkdir -p "$out_dir"

  out_path="${out_dir}/evaluation_result_100_10000_transformer_${idx}.npz"

  CUDA_VISIBLE_DEVICES=0 python /path/to/transformer_benchmark.py \
    --test_data "$i" \
    --device cuda \
    --weights "/path/to/model.safetensors" \
    --save_path "$out_path"
done