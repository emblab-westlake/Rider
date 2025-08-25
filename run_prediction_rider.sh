
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
#RIDER_PATH=$(dirname ${SCRIPT_PATH})

cd $SCRIPT_PATH

# source /opt/miniforge3/bin/activate rider
INPUT_DIR=/root/gaoyang/westlake_emblab/Rider/test_data  #test_data
output_path=/root/gaoyang/westlake_emblab/Rider/test_data/test_results #results
weights=$SCRIPT_PATH/checkpoint/checkpoint-102000/model.safetensors
Rider_pdb_database=/usr/commondata/public/gaoyang/Rider_pdb_database/database

# Print for debugging
echo "Using weights from: $weights"

for i in $INPUT_DIR/*faa
do

base=$(basename ${i})
File_path=$(dirname ${i})
out_dir=${output_path}/${base}
mkdir -p ${out_dir}

CUDA_VISIBLE_DEVICES=6 \
    python predict_pipline.py \
        -i ${i} \
        -t 32 \
        -w ${weights} \
        -b 64 \
        -o ${output_path} \
        --predict_structure \
        --sequence_length 1024 \
        --structure_align_enabled \
        --rdrp_structure_database ${Rider_pdb_database} \
        --prob_threshold 50 \
        --top_n_mean_prob 1 \
        --alignment-type 1

done
