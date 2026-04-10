#!/bin/bash
set -e

PYTHON=python

ESM_PATH=esm2_t12_35M_UR50D
PRETRAIN_PATH=$ESM_PATH/model.safetensors
OUTPUT_DIR=training_data
DEVICE=cuda:1
SEED=42

TOKENIZE_SCRIPT=step1_tokenize_fasta_to_pt.py
BUILD_SCRIPT=step2_build_train_embeddings.py

RAW_DIR=train_test_set

VIRUS_FASTA=$RAW_DIR/train_data.fasta
NON_VIRUS_FASTA=$RAW_DIR/input_ids_nonvirus.fasta

PROTEASE_FASTA=$RAW_DIR/other_rna_protease_train_data.fasta
CAPSID_FASTA=$RAW_DIR/other_rna_capsid_train_data.fasta
HELICASE_FASTA=$RAW_DIR/other_rna_helicase_train_data.fasta

VIRUS_PT=$RAW_DIR/train_data.pt
NON_VIRUS_PT=$RAW_DIR/input_ids_nonvirus.pt
PROTEASE_PT=$RAW_DIR/other_rna_protease_train_data.pt
CAPSID_PT=$RAW_DIR/other_rna_capsid_train_data.pt
HELICASE_PT=$RAW_DIR/other_rna_helicase_train_data.pt

# 1) tokenization
$PYTHON $TOKENIZE_SCRIPT \
  --esm_path $ESM_PATH \
  --input_fasta $VIRUS_FASTA \
  --output_pt $VIRUS_PT

$PYTHON $TOKENIZE_SCRIPT \
  --esm_path $ESM_PATH \
  --input_fasta $NON_VIRUS_FASTA \
  --output_pt $NON_VIRUS_PT

$PYTHON $TOKENIZE_SCRIPT \
  --esm_path $ESM_PATH \
  --input_fasta $PROTEASE_FASTA \
  --output_pt $PROTEASE_PT

$PYTHON $TOKENIZE_SCRIPT \
  --esm_path $ESM_PATH \
  --input_fasta $CAPSID_FASTA \
  --output_pt $CAPSID_PT

$PYTHON $TOKENIZE_SCRIPT \
  --esm_path $ESM_PATH \
  --input_fasta $HELICASE_FASTA \
  --output_pt $HELICASE_PT

# 2) normal
$PYTHON $BUILD_SCRIPT \
  --mode normal \
  --device $DEVICE \
  --esm_path $ESM_PATH \
  --pretrain_path $PRETRAIN_PATH \
  --virus_input_ids_path $VIRUS_PT \
  --non_virus_input_ids_path $NON_VIRUS_PT \
  --output_dir $OUTPUT_DIR \
  --output_prefix train_normal \
  --neg_pos_ratio 5 \
  --seed $SEED

# 3) hard negative
$PYTHON $BUILD_SCRIPT \
  --mode hard_negative \
  --device $DEVICE \
  --esm_path $ESM_PATH \
  --pretrain_path $PRETRAIN_PATH \
  --virus_input_ids_path $VIRUS_PT \
  --non_virus_input_ids_path $NON_VIRUS_PT \
  --extra_neg_protease $PROTEASE_PT \
  --extra_neg_capsid $CAPSID_PT \
  --extra_neg_helicase $HELICASE_PT \
  --output_dir $OUTPUT_DIR \
  --output_prefix train_hardneg_3extra \
  --neg_pos_ratio 5 \
  --seed $SEED