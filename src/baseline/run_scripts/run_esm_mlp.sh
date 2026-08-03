#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:" >&2
  echo "  $0 EMBEDDINGS_DIR OUTPUT_DIR WEIGHTS" >&2
  echo "  $0 EMBEDDINGS_DIR OUTPUT_DIR --untrained-control" >&2
}

if [[ "$#" -ne 3 ]]; then
  usage
  exit 2
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASELINE_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
EMBEDDINGS_DIR="$1"
OUT_ROOT="$2"
MODE_OR_WEIGHTS="$3"
PYTHON="${PYTHON:-python}"
DEVICE="${DEVICE:-cuda}"
SEED="${SEED:-42}"

if [[ ! -d "$EMBEDDINGS_DIR" ]]; then
  echo "Embedding directory not found: $EMBEDDINGS_DIR" >&2
  exit 1
fi

if [[ "$MODE_OR_WEIGHTS" == "--untrained-control" ]]; then
  MODEL_LABEL="ESM2_untrained"
  MODEL_ARGS=(--untrained-control --seed "$SEED")
elif [[ -f "$MODE_OR_WEIGHTS" ]]; then
  MODEL_LABEL="ESM_mlp"
  MODEL_ARGS=(--weights "$MODE_OR_WEIGHTS")
else
  echo "Checkpoint not found: $MODE_OR_WEIGHTS" >&2
  exit 1
fi

mkdir -p "$OUT_ROOT"
shopt -s nullglob
embedding_files=("$EMBEDDINGS_DIR"/*.pt)
if [[ "${#embedding_files[@]}" -eq 0 ]]; then
  echo "No .pt embedding files found in: $EMBEDDINGS_DIR" >&2
  exit 1
fi

for i in "${embedding_files[@]}"; do
  echo "Processing file: $i"

  base="$(basename "$i")"
  stem="${base%.pt}"

  # Extract the set index from the filename, e.g. set1 -> 1
  idx="$(echo "$stem" | sed -n 's/.*_set\([0-9]\+\)_.*/\1/p')"
  idx="${idx:-X}"

  # Create output directory for each set
  out_dir="${OUT_ROOT}/set${idx}"
  mkdir -p "$out_dir"

  out_path="${out_dir}/evaluation_result_100_10000_${MODEL_LABEL}_${idx}.npz"

  "$PYTHON" "$BASELINE_DIR/esm_mlp_benchmark.py" \
    --embeddings_pt "$i" \
    "${MODEL_ARGS[@]}" \
    --device "$DEVICE" \
    --save_path "$out_path"
done
