# 🧬 Rider

**Rider** enables fast identification of known and novel RNA viruses from large volumes of metatranscriptomic sequencing data.  
It integrates sequence classification, structure predictcion, and structure alignment into a streamlined pipeline.

We are keeping updating this project...


## 🏗️ Architecture

![Rider Workflow](docs/architecture/Rider_workflow.png)
*Figure 1: Rider pipeline workflow showing the complete data processing flow from input sequences to final virus identification.*

## 🚀 Installation
Two supported ways to prepare the Rider environment:

**A. Developer / reproducible setup (recommended for development)**

This keeps your original env.yaml-based workflow. It recreates the exact conda environment used during development and installs Python packages from requirements.txt.

You can use `git clone` and `conda` to set up the environment:

```bash
# Clone the repository
git clone https://github.com/emblab-westlake/Rider.git
cd Rider

# Create conda environment
conda env create -f env.yaml
conda activate rider

# Install Python-only dependencies (pip)
pip install -r requirements.txt

# Install Rider package in editable mode for development
pip install -e .
```
Notes:

- `env.yaml` installs system/bio tools (mmseqs2, diamond, blast, hmmer, prodigal, entrez-direct) from conda-forge / bioconda so they are available on PATH.
- `requirements.txt` lists many heavy GPU and system packages (torch, triton, deepspeed). Installing them via pip in some systems may fail or produce suboptimal GPU behavior — see the alternative installation below for a safer approach.

**B. Manual Conda Setup (GPU/CPU control; safer for end users; preferred for running pipelines)**

Use conda to install PyTorch (choose CPU or GPU build) and binary tools, then install the Rider package with pip. This avoids pip attempting to build or incorrectly install GPU binaries.
CPU-only example:
```bash
conda create -n rider python=3.10 -y
conda activate rider

# Install PyTorch (CPU-only) via official PyTorch channel
conda install -c pytorch pytorch cpuonly -y

# Install bioinformatics binaries via bioconda/conda-forge
conda install -c conda-forge -c bioconda mmseqs2 foldseek diamond blast hmmer prodigal entrez-direct -y

# Install other Python dependencies (pip) and Rider (editable)
pip install -r requirements.txt --no-deps
pip install -e .
```
GPU example (adjust cudatoolkit to your GPU/cuda compatibility):
```bash
conda create -n rider-gpu python=3.10 -y
conda activate rider-gpu

# Install PyTorch with CUDA support (example uses cudatoolkit=11.8)
conda install -c pytorch pytorch torchvision torchaudio cudatoolkit=11.8 -y

# Install bioinformatics binaries
conda install -c conda-forge -c bioconda mmseqs2 foldseek diamond blast hmmer prodigal entrez-direct -y

# Install other Python deps without re-installing torch/triton/deepspeed (use --no-deps)
pip install -r requirements.txt --no-deps
pip install -e .
```

**Quick verification**

After installation, ensure binaries are available:
```bash
which mmseqs
mmseqs -h
```
Then test Rider CLI (after pip install -e .):
```bash
rider-predict -h
```

## 📁  Step 2. Download and prepare the Rider Dependency
### 📦 Dependency
Rider depends on the following:

- ESM2 (sequence embeddings)
- ESMFold (structure prediction)
- Foldseek (structure alignment)
- Rider structure database (RNA viral protein strucutre reference) (comming soon...)

This database is required for structural alignment using Foldseek (step 6/7).

📥 How to prepare:
1. Download the prebuilt dependencies and database manually (Zenodo: https://doi.org/10.5281/zenodo.15742756) or follow your internal distribution process. We are preparing this...
2. Place it under the Rider folder like this:

✅ Required layout:
```sh
Rider/
└── Submodule/
    └── esmfold_v1/    # esmfold
    └── esm2_t12_35M_UR50D # ESM
    └── foldseek #foldseek binaries (e.g., foldseek-linux-*)
└── Rider_pdb_database/
    └── database/    # <- Foldseek database files go here
```

Alternatively, pass a custom path using:
```sh
--rdrp_structure_database /your/custom/path/database
```
Note: Foldseek and mmseqs2 are external binaries. Rider attempts to add submodule/foldseek/bin to PATH, but mmseqs2 usually needs separate installation (e.g., conda: `conda install -c bioconda mmseqs2`). Verify with `which mmseqs` and `mmseqs --version`.

## 🧪 Quick run example
You can run the pipeline using the provided shell script `Rider/run_rider_predict.sh`:

```sh
# Activate the conda environment
# source /root/miniconda3/bin/activate rider
SCRIPT_PATH=$(dirname $(readlink -f "$0"))
# Set input and output paths
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
        --device cpu \ #or --device gpu
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
```
Arguments explained

- `-i`, `--input_faa` (str, required)
Path to input FASTA file. Each record should be one protein sequence.

- `-w`, `--weights` (str, required)
Path to the classification model weights (safetensors). Default in code: checkpoint/checkpoint-102000/model.safetensors.

- `-b`, `--batch_sizes` (int, default=64)
Batch size for tokenization / feature extraction. Adjust based on GPU memory (≤ 64 suggested for <16GB GPU).

- `-t`, `--threads` (int, default=16)
Number of CPU threads.

- `-o`, `--output_dir` (str, default=./)
Output root directory. The pipeline creates a per-input subdirectory with intermediate and final results.

- `--sequence_length` (int, default=1024)
Tokenizer maximum length. Increase only if necessary and if model supports it.

- `--predict_structure` (flag)
Enable structure prediction (ESMFold). If not set, structure prediction steps are skipped.

- `--structure_model_path` (str)
Path to a custom structure model. If not provided, submodule/esmfold_v1 is used.

- `--threshold` (float, default=0.8)
Classification score threshold to call positive.

- `--device` (str, default="cuda")
Computation device (e.g., "cuda", "cuda:0", "cpu").

- `--negative_sample_path` (str)
Path to negative sample FASTA used to pad batches. Default: databases/false256.faa.

- `--structure_align_enabled` (flag)
Enable Foldseek structural alignment (requires Foldseek binary and Rider PDB database).

- `--rdrp_structure_database` (str)
Path to Foldseek database directory. Required if --structure_align_enabled is set.

- `--alignment-typ`e (int, default=1)
Foldseek alignment type parameter (passed to the alignment runner).

- `--prob_threshold` (int, default=50)
Homology probability threshold (percentage) used to filter Foldseek results.

- `--top_n_mean_prob` (int, default=1)
Number of top hits to average when computing homology probability. Higher values make validation stricter.

## 🔗Cite us


<!-- ## 🔗 Merge Foldseek top-hits with taxonomy

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

- --m8_file : Path to the Foldseek .m8 file (tab-separated, columns: query, target, identity, aln_len, matches, mismatches, q_start, q_end, t_start, t_end, evalue, tm_score).
- --taxo_file : Path to a taxonomy TSV file that must contain an Accession column.
- --output_file : Path to write the merged output (TSV).

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
``` -->

