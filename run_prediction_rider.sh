
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
#RIDER_PATH=$(dirname ${SCRIPT_PATH})

cd $SCRIPT_PATH

source /root/miniconda3/bin/activate LRVM_gy
INPUT_DIR=/home/gaoyang/dataset/yuanlin/find_new_vir/complete_proteins  #yuanlin_sewage
output_path=/usr/commondata/public/gaoyang/software/rider/RNA_virus_project/yuanlinAS #yuanlin_AS
weights=$SCRIPT_PATH/checkpoint/checkpoint102000/model.safetensors

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
        -w ${weights} \
        -b 256 \
        -o ${output_path} \
        --predict_structure \
        --sequence_length 1024 \
        --structure_align_enabled \
        --rdrp_structure_database "/usr/commondata/public/gaoyang/Rider_pdb_database/database" \
        --prob_threshold 50 \
        --top_n_mean_prob 1 \
        --alignment-type 1

done
