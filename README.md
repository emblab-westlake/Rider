# 🧬 Rider

**Rider** enables fast identification of known and novel RNA viruses from large volumes of metatranscriptomic sequencing data.  
It integrates sequence classification, structure predictcion, and structure alignment into a streamlined pipeline.


## 🏗️ Architecture

![Rider Workflow](docs/architecture/Rider_workflow_3.png)
*Figure 1: Rider pipeline workflow showing the complete data processing flow from input sequences to final virus identification.*

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

## 📁  Step 2. Download and prepare the Rider RDRP structure database

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

## 🔗 Merge Foldseek top-hits with taxonomy

A helper script is provided to merge Foldseek m8 results (one-line-per-alignment format)
with a taxonomy table. The script selects the top hit per query according to the `tm_score`
and attaches taxonomy fields (matched by accession, with optional version suffix removed).

Script path: `Rider/utils/rider_merge_taxonomy.py`

Usage:
```bash
python /root/gaoyang/westlake_emblab/Rider/utils/rider_merge_taxonomy.py \
    --m8_file /path/to/foldseek_results.m8 \
    --taxo_file /path/to/taxonomy.tsv \
    --output_file /path/to/output_top_hit_taxonomy.tsv
```

Arguments:

- --m8_file : Path to the Foldseek .m8 file (tab-separated, columns: query, target, identity, aln_len, matches, mismatches, q_start, q_end, t_start, t_end, evalue, tm_score). Default in script:
/usr/commondata/public/gaoyang/foldseek_search_dir/rider_gt200_all_mapping_aln.m8
- --taxo_file : Path to a taxonomy TSV file that must contain an Accession column. Default in script:
/home/gaoyang/Rider/utils/taxo_Rider_train_dat_parsed_taxonomy.tsv
- --output_file : Path to write the merged output (TSV). Default in script:
/usr/commondata/public/gaoyang/foldseek_search_dir/foldseek_top_hit_taxonomy_lucaprot_gt200_all_mapping_aln.tsv

What the script does:

1. Reads the m8 file and converts the tm_score column to numeric.
2. Picks the top hit per query by descending tm_score.
3. Removes version suffix from target (e.g. .1) before matching to accession.
4. Reads the taxonomy file, expects an Accession column, and deduplicates by Accession.
5. Left-joins the top hits with taxonomy fields using cleaned target -> Accession.
6. Outputs a TSV with columns: query, target, plus all taxonomy columns except Accession.

Notes and tips:

- The script will exit with an error if either input file is missing or if the taxonomy file lacks an Accession column.
- If your m8 file or taxonomy file is large, ensure sufficient memory is available. Consider pre-filtering if necessary.
- The script strips a trailing dot followed by digits from the target accession (regex: \.\d+$). Modify the regex in the script if your accession format differs.
- To make the script executable and run it directly:
```bash
chmod +x /Rider/utils/rider_merge_taxonomy.py
/root/gaoyang/westlake_emblab/Rider/utils/rider_merge_taxonomy.py --m8_file ... --taxo_file ... --output_file ...
```
Example:
```bash
python Rider/utils/rider_merge_taxonomy.py \
--m8_file /usr/commondata/public/gaoyang/foldseek_search_dir/rider_gt200_all_mapping_aln.m8 \
--taxo_file Rider/utils/taxo_Rider_train_dat_parsed_taxonomy.tsv \
--output_file /foldseek_search_dir/foldseek_top_hit_taxonomy_lucaprot_gt200_all_mapping_aln.tsv
```
