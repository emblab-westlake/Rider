
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
#RIDER_PATH=$(dirname ${SCRIPT_PATH})

cd $SCRIPT_PATH

# source /opt/miniforge3/bin/activate rider
INPUT_DIR=$SCRIPT_PATH/test_data  #test_data
output_path=$SCRIPT_PATH/test_data/test_results #results
weights=$SCRIPT_PATH/checkpoint/checkpoint-44000/model.safetensors
Rider_pdb_database=$SCRIPT_PATH/Rider_pdb_database/database

# Print for debugging
echo "Using weights from: $weights"

for i in $INPUT_DIR/*faa
do

base=$(basename ${i})
File_path=$(dirname ${i})
out_dir=${output_path}/${base}
mkdir -p ${out_dir}

CUDA_VISIBLE_DEVICES=0 \
    python predict_pipline.py \
        -i ${i} \
        -t 32 \
        -w ${weights} \
        -b 64 \
        -o ${output_path} \
        --sequence_length 1024 \
        --predict_structure \
        --structure_align_enabled \
        --rdrp_structure_database ${Rider_pdb_database} \
        --prob_threshold 50 \
        --top_n_mean_prob 2 \
        --alignment-type 1

done
