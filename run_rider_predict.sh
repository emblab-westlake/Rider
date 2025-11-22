#!/bin/bash
#conda activate rider

# get the directory of the current script
SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "${SCRIPT}")

# switch to the script directory to ensure relative paths are correct
cd "$SCRIPT_PATH"

# input, output paths and model weights
INPUT_DIR="$SCRIPT_PATH/test_data"                          # input dir
OUTPUT_PATH="$SCRIPT_PATH/test_data/test_results"           # output dir
WEIGHTS="$SCRIPT_PATH/checkpoint/checkpoint-44000/model.safetensors"
RDRP_DB="/usr/commondata/public/gaoyang/Rider_pdb_database/Rider_pdb_database/database" #change this to your own path
SUBMODULE_DIR="$SCRIPT_PATH/submodule"

# Debug output
echo "Using weights from: $WEIGHTS"
echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_PATH"

# loop through all .faa files in the input directory
for i in "$INPUT_DIR"/*faa; do
    base=$(basename "$i")
    File_path=$(dirname "$i")
    out_dir="${OUTPUT_PATH}/${base}"
    mkdir -p "${out_dir}"

    echo "Processing: $i"

    # run rider-predict 
    # Set CUDA_VISIBLE_DEVICES to specify which GPU to use
    CUDA_VISIBLE_DEVICES=4 \
    rider-predict \
        -i "$i" \
        -t 32 \
        -w "$WEIGHTS" \
        -b 64 \
        --device cpu \
        -o "$OUTPUT_PATH" \
        --submodule_dir "$SUBMODULE_DIR" \
        --predict_structure \
        --sequence_length 1024 \
        --structure_align_enabled \
        --rdrp_structure_database "$RDRP_DB" \
        --prob_threshold 50 \
        --top_n_mean_prob 2 \
        --alignment-type 1

    echo "Finished: $i"
    echo "--------------------------------------------"
done