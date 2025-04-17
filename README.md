# Rider
Rider enables fast identification of known and novel RNA rirus from large volum of metatrasncriptomic sequencing data.

## Installation
You can use `git clone` or use `conda install` to set up environments.
```sh
#install via git clone
git clone https://github.com/emblab-westlake/Rider.git
cd Rider
conda env create -f environment.yaml
cd submodel
tar xvzf foldseek-linux-avx2.tar.gz
```
Rider need dependencies include `ESM2_t12` and `ESMFold_v1` model. Make sure these two directories are under the Rider fold.

## Easy example
Your can simply run `run_prediction_rider.sh`:
```sh
conda activate rider
#recommend to use single gpu for inference
INPUT_DIR= #set up your gene .faa files
CUDA_VISIBLE_DEVICES=0 \
for i in $INPUT_DIR/*INF*faa
do

base=$(basename ${i})
File_path=$(dirname ${i})
mkdir -p ${out_dir}

CUDA_VISIBLE_DEVICES=4 \
    python predict_pipline.py \
        -i ${i} \
        -w ${weights} \
        -b 256 \
        -o ${output_path} \
        --predict_structure \
        --sequence_length 1024 \
        --structure_align_enabled \
        --rdrp_structure_database "/usr/commondata/public/gaoyang/Rider_pdb_database/database"

done
```

## Testing on the benchmark
Simply run the command `bash run_testing.sh`. It is recommended to use single GPU for inference.
Your can also change parameter in the file `run_finetune.sh`:
```sh

CUDA_VISIBLE_DEVICES=4 \
python benchmark_test.py \
-w <abs path to the finetuned saved safe.tensor> \
-b <batch size>
```

## How to Fine-tune
In file finetune.py, we can change pretrain weights path, positive data path and negative data path.

```
The postigve data path: 
```sh
virus_input_ids_path = "/usr/commondata/public/gaoyang/dataset/LRVM_train_test_set/filterd_rdrp_positive_finetune_19634_diamond_mmseq.pt"
```
The negative data path:
```sh
#We recommend to use uniref50 after filtering 
non_virus_input_ids_path = "/usr/commondata/public/gaoyang/dataset/LRVM_train_test_set/filterd_non_rdrp_negative_finetune_50w.pt"
```

## How to change ratio of positive and negative
In file `model.py`, you can change weights when calculating loss \
`self.celoss = torch.nn.CrossEntropyLoss(ignore_index=-100, weight=torch.tensor([1.0,10.0]))` \
Here we use 1:10, where label0: label1 = 1:10, as the size of label 1 data are relatively smaller.