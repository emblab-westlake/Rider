#!/bin/bash

# Run Rider benchmark on multiple embedding files
cd /path/to/baseline/scripts

# Output root directory for results
out_root="path/to/save/results"

# Iterate over embedding files
for i in /path/to/test_set/ESM2_35M_embeddings/*.pt; do
  echo "Processing file: $i"

  base="$(basename "$i")"
  stem="${base%.pt}"

  idx="$(echo "$stem" | sed -n 's/.*_set\([0-9]\+\)_.*/\1/p')"
  idx="${idx:-X}"

  out_dir="${out_root}/set${idx}"
  mkdir -p "$out_dir"

  out_path="${out_dir}/evaluation_result_rider_${idx}.npz"

  CUDA_VISIBLE_DEVICES=0 python /path/to/rider_benchmark.py \
    --embeddings_pt "$i" \
    --weights "/path/to/model.safetensors" \
    --device cuda \
    --save_path "$out_path"
done