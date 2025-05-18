# 🧬 Rider

**Rider** enables fast identification of known and novel RNA viruses from large volumes of metatranscriptomic sequencing data.  
It integrates sequence classification, structure prediction, and structure alignment into a streamlined pipeline.

---

## 🚀 Installation
You can use `git clone` and `conda` to set up the environment:

```bash
# Clone the repository
git clone https://github.com/emblab-westlake/Rider.git
cd Rider

# Create conda environment
conda env create -f environment.yaml

# Unpack Foldseek
cd submodule
tar xvzf foldseek-linux-avx2.tar.gz

# Unpack ESM2 model (if compressed)
cd ..
gunzip esm2_t12_35M_UR50D.gz
```
✅ Make sure the following two folders exist under the Rider root:

- `esm2_t12_35M_UR50D/`
- `esmfold_v1/`

## 📦 Dependency
Rider depends on the following:

- ESM2 (sequence embeddings)
- ESMFold (structure prediction)
- Foldseek (structure alignment)

All Python dependencies are included in environment.yaml.

## 📁  Step 2. Download and prepare the Foldseek RDRP structure database

This database is required for structural alignment using Foldseek (step 6/7).

📥 How to prepare:
1. Download the Foldseek-format database manually (link to be provided).
2. Place it under the Rider folder like this:
```sh
Rider/
└── Rider_pdb_database/
    └── database/    # <- Foldseek database files go here
```

Alternatively, pass a custom path using:
```sh
--rdrp_structure_database /your/custom/path/database
```


## 🧪 Easy Example
You can run the pipeline using the provided shell script `run_prediction_rider.sh`:

```sh
# Activate the conda environment
source /root/miniconda3/bin/activate LRVM_gy

# Set input and output paths
INPUT_DIR=/home/gaoyang/dataset/yuanlin/find_new_vir/complete_proteins
OUTPUT_DIR=/usr/commondata/public/gaoyang/software/rider/RNA_virus_project/yuanlinAS
SCRIPT_PATH=$(dirname $(readlink -f "$0"))
WEIGHTS=$SCRIPT_PATH/checkpoint/checkpoint102000/model.safetensors
```sh

```sh
for i in $INPUT_DIR/*faa
do
    base=$(basename ${i})
    out_dir=${OUTPUT_DIR}/${base}
    mkdir -p ${out_dir}

    CUDA_VISIBLE_DEVICES=6 \
    python predict_pipline.py \
        -i ${i} \
        -w ${WEIGHTS} \
        -b 256 \
        -o ${OUTPUT_DIR} \
        --predict_structure \
        --sequence_length 1024 \
        --structure_align_enabled \
        --rdrp_structure_database "$SCRIPT_PATH/Rider_pdb_database/database" \
        --prob_threshold 50 \
        --top_n_mean_prob 1 \
        --alignment-type 1
done

```