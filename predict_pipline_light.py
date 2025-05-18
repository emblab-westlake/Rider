"""
------------------------------------------
Script Name: predict_pipline.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - The main workflow of Rider
Usage: 
    - python predict_pipline.py -h
Dependencies: 
    - Python 3.9 or higher
    - pandas >= 1.3.0
    - numpy >= 1.21.0
Notes:
    - Enjoy yourself!
    - Output directory will be created if it does not exist
------------------------------------------
"""

import os
import time
import json
import psutil
import torch
import numpy as np
from torch import nn
from torch.utils.data import DataLoader
from transformers import AutoTokenizer, EsmModel
from transformers.utils import ModelOutput
from Bio import SeqIO
from safetensors.torch import load_model
from datasets import Dataset
import argparse
import subprocess
import logging
from src.predict_structure_package import load_structure_model, process_faa_files
from src.mmcluster import mmseqs2_clustering  
from src.classification import classify_genes
from src.feature_extraction import extract_features
from src.structure_align import Structure_aligned
from src.filter_prob_taxonomy import process_fixed_prob_out_folders
from src.parse_result import process_final_result  

# logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# parse arg
def parse_args():
    parser = argparse.ArgumentParser(description="Protein sequence classification and structure prediction pipeline.")
    parser.add_argument("-i", "--input_faa", required=True, type=str, 
                        help="Path to the input FASTA file.")
    parser.add_argument("-w", "--weights", required=True, type=str, 
                        help="Path to the model weights.")
    parser.add_argument("-b", "--batch_sizes", type=int, default=64, 
                        help="Batch size for processing.")
    parser.add_argument("-t", "--threads", type=int, default=16, 
                        help="Number of CPU threads to use.")
    parser.add_argument("-o", "--output_dir", type=str, default="./", 
                        help="Output directory.")
    parser.add_argument("--sequence_length", type=int, default=1024, 
                        help="Protein sequence truncation length.")
    parser.add_argument("--predict_structure", action="store_true", 
                        help="Enable structure prediction.")
    parser.add_argument("--structure_model_path", type=str, 
                        help="Path to the structure prediction model.")
    parser.add_argument("--threshold", type=float, default=0.5, 
                        help="Threshold value for classification (default: 0.5).")
    parser.add_argument("--device", type=str, default="cuda", 
                        help="Device to use for computation (default: 'cuda').")
    parser.add_argument("--negative_sample_path", type=str, default="/home/gaoyang/dataset/test_pipeline/test_data/false256.faa", 
                        help="Path to the negative samples file.")
    parser.add_argument("--structure_align_enabled", action="store_true", 
                        help="Enable Foldseek alignment task.")
    parser.add_argument("--rdrp_structure_database", type=str, required=False, 
                        help="Path to Foldseek database.")
    parser.add_argument("--alignment-type", type=int, default=1, 
                        help="Foldseek alignment type parameter.")
    parser.add_argument("--prob_threshold", type=int, default=50, 
                        help="Homologous probability threshold default 50.")
    parser.add_argument("--top_n_mean_prob", type=int, default=1, 
                        help="Usually the at least top 1 query above the homo prob threshold, you can choose 2 or more to obtain more strict results.")
    return parser.parse_args()

# model
class LRVMForMaskedLM(nn.Module):
    def __init__(self, model_path):
        super(LRVMForMaskedLM, self).__init__()
        self.esm = EsmModel.from_pretrained(model_path)
        self.ln = nn.Linear(480, 33)
        self.celoss = nn.CrossEntropyLoss(ignore_index=-100)
    
    def forward(self, input_ids):
        x = self.esm(**{'input_ids': input_ids})['last_hidden_state']
        logits = x
        loss = torch.tensor(0.0) 
        return ModelOutput({'logits': logits, 'loss': loss})

class LRVMForClf(nn.Module):
    def __init__(self, esm_model_path, pretrain_path, clf_class):
        super(LRVMForClf, self).__init__()
        self.LRVM = LRVMForMaskedLM(esm_model_path)
        load_model(self.LRVM, pretrain_path, strict=False)
        self.ln = nn.Linear(480, clf_class)
        self.celoss = nn.CrossEntropyLoss(ignore_index=-100)
    
    def forward(self, input_ids):
        x = self.LRVM(**{'input_ids': input_ids})['logits'][:, 0, :]
        features = x
        x = self.ln(features)
        logits = x
        return ModelOutput({'logits': logits, 'features': features})

class SimpleClassifier(nn.Module):
    def __init__(self, input_dim=480, output_dim=2, hidden_dim=1420, nhead=8, num_layers=1, dropout=0.3):
        super(SimpleClassifier, self).__init__()
        encoder_layer = nn.TransformerEncoderLayer(d_model=input_dim, nhead=nhead, dim_feedforward=hidden_dim, dropout=dropout)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.linear = nn.Linear(input_dim, output_dim)
    
    def forward(self, input_ids=None, labels=None):
        transformer_output = self.transformer_encoder(input_ids)
        if transformer_output.dim() == 3:
            first_token_output = transformer_output[:, 0, :]
        elif transformer_output.dim() == 2:
            first_token_output = transformer_output
        else:
            raise ValueError(f"Unexpected transformer output dimensions: {transformer_output.dim()}")
        logits = self.linear(first_token_output)
        loss = None
        if labels is not None:
            loss_fct = nn.CrossEntropyLoss()
            loss = loss_fct(logits, labels)
        return ModelOutput({'loss': loss, 'logits': logits, 'features': first_token_output})

def load_and_tokenize_data(input_faa, tokenizer, sample_source_path, max_length=1024, batch_size=256):
    """
    Load input samples and tokenize them. If the number of input sequences is less than batch_size,
    pad with sequences from a sample source.

    Args:
        input_faa (str): Path to the input FASTA file containing sample sequences.
        tokenizer (Callable): Function or tokenizer to tokenize the sequences.
        sample_source_path (str): Path to the sample source FASTA file for padding.
        max_length (int): Maximum sequence length for truncation (default: 1024).
        batch_size (int): Target batch size. Padding is applied if input count is less.

    Returns:
        dict: A dictionary with the following keys:
            - "encoded_sequences" (Dataset): Tokenized and padded sequence dataset.
            - "padded_indices" (list): Indices of the padded sequences in the dataset.
        list: Original sequence dictionary from input_faa, with IDs as keys and sequences as values.
    """
    # Load input sequences
    seq_store = {}
    for record in SeqIO.parse(input_faa, 'fasta'):
        seq_store[str(record.id)] = str(record.seq)
    data_lst = list(seq_store.values())
    logging.info(f"Number of input sequences: {len(data_lst)}")

    # Load sample source sequences
    sample_source = []
    for record in SeqIO.parse(sample_source_path, 'fasta'):
        sample_source.append(str(record.seq))
    logging.info(f"Number of sequences in sample source: {len(sample_source)}")

    # Check if padding is needed
    padded_indices = []
    num_inputs = len(data_lst)
    if num_inputs < batch_size:  # Only pad if inputs are fewer than batch size
        padding_size = batch_size - num_inputs
        if len(sample_source) >= padding_size:
            # Select padding sequences from the sample source
            padded_samples = sample_source[:padding_size]
            data_lst.extend(padded_samples)
            padded_indices = list(range(num_inputs, num_inputs + len(padded_samples)))
            logging.info(f"Padded {len(padded_samples)} samples to match batch size.")
        else:
            logging.warning(f"Not enough samples in the source to pad. Expected {padding_size}, found {len(sample_source)}.")

    # Create a HuggingFace Dataset
    ds = Dataset.from_dict({'data': data_lst})

    # Tokenization
    encoded_sequences = ds.map(
        lambda examples: tokenizer(
            examples['data'], truncation=True, max_length=max_length, padding="max_length", return_tensors="pt"
        ),
        batched=True,
        remove_columns='data'
    )
    encoded_sequences.set_format('pt')  # Set format to PyTorch tensors

    return {
        "encoded_sequences": encoded_sequences,
        "padded_indices": padded_indices
    }, seq_store


def save_predicted_genes(output_dir, file_name, seq_store_pos, target_positive_record_id, target_negative_record_id, input_faa):
    """
    Save predicted gene results to the specified directory, including:
    - All prediction results file
    - FASTA file for RNA virus candidate sequences
    - FASTA file for non-RNA virus sequences

    Args:
    - output_dir: Output directory
    - file_name: Base name for output files (derived from input file)
    - seq_store_pos: Dictionary of sequences {id: sequence}
    - target_positive_record_id: List of sequence IDs predicted as positive
    - target_negative_record_id: List of sequence IDs predicted as negative
    - input_faa: Path to the original input FASTA file
    """
    start_time = time.time()

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save all prediction results to a text file
    results_file = os.path.join(output_dir, f"{file_name}_Rider_predicted_results.txt")
    with open(results_file, "w") as f:
        for i in target_positive_record_id:
            f.write(i + "\t" + seq_store_pos[i] + "\t" + str(1) + "\n")
        for j in target_negative_record_id:
            f.write(j + "\t" + seq_store_pos[j] + "\t" + str(0) + "\n")
    logging.info(f"Predicted results saved to: {results_file}")

    # Save RNA virus candidate sequences in FASTA format
    positive_fasta_file = os.path.join(output_dir, f"{file_name}_Rider_predicted_RNA_Virus_potential_candidates.faa")
    with open(positive_fasta_file, "w") as f:
        for i in target_positive_record_id:
            f.write(">Rider_" + i + "\n")
            f.write(seq_store_pos[i] + "\n")
    logging.info(f"RNA Virus candidates saved to: {positive_fasta_file}")

    # Save non-RNA virus sequences in FASTA format
    negative_fasta_file = os.path.join(output_dir, f"{file_name}_Rider_predicted_nonRNA.faa")
    with open(negative_fasta_file, "w") as f:
        for j in target_negative_record_id:
            f.write(">" + j + "\n")
            f.write(seq_store_pos[j] + "\n")
    logging.info(f"Non-RNA sequences saved to: {negative_fasta_file}")

    # Log completion info
    logging.info(f"{input_faa} done.")
    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Saving predicted genes completed in {elapsed_time:.2f} seconds.")

def main_light():
    overall_start_time = time.time()
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss / 1024 / 1024  # MB

    # ======= Up to Step 3 =======

    # Get the path of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set the absolute path to foldseek_dir
    foldseek_dir = os.path.join(script_dir, "submodule", "foldseek", "bin")
    
    # Dynamically update the PATH environment variable
    os.environ["PATH"] = f"{foldseek_dir}:{os.environ.get('PATH', '')}"

    # Optional: print updated PATH to verify
    print(f"Updated PATH: {os.environ['PATH']}")
    
    args = parse_args()
    start_time = time.time()
    overall_start_time = start_time  # Record the start time of the entire pipeline
    
    # Path configurations
    file_name = os.path.basename(args.input_faa)
    tmp_dir = os.path.join(args.output_dir, file_name, f"{file_name}_intermediate")
    tmp_tensor_dir = os.path.join(tmp_dir, "tmp_tensors")
    embeddings_tensor_path = os.path.join(tmp_tensor_dir, f"{file_name}_features_embeddings.pt")
    tmp_mmseq_dir = os.path.join(tmp_dir, "tmp_mmseqs")
    tmp_mmseq_tmp_dir = os.path.join(tmp_mmseq_dir, "tmp")
    predicted_rdrp_fasta_path = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_RNA_Virus_potential_candidates.faa")
    
    # Model and resource paths
    esmfold_dir = os.path.join(script_dir, "esmfold_v1")
    esmt12_dir = os.path.join(script_dir, "esm2_t12_35M_UR50D")
    known_rdrp = os.path.join(script_dir, "data", "NCBI_RNA_virus_refseq.fasta")

    # Output file paths
    results_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_results.txt")
    positive_fasta_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_RNA_Virus_potential_candidates.faa")
    negative_fasta_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_nonRNA.faa")
    
    # Ensure required directories exist
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(tmp_tensor_dir, exist_ok=True)
    os.makedirs(tmp_mmseq_dir, exist_ok=True)
    os.makedirs(tmp_mmseq_tmp_dir, exist_ok=True)
    
    # Set up log file path
    log_file_path = os.path.join(tmp_dir, "rider_pipeline.log")
    os.makedirs(args.output_dir, exist_ok=True)  # Ensure the output directory exists
    
    # Configure logging to both file and console
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Clear existing handlers to avoid duplicates
    if logger.hasHandlers():
        logger.handlers.clear()

    # Define log formatter
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # File handler: write logs to file
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler: output logs to console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Verify logging
    logging.info(f"Logging to file: {log_file_path}")
    logging.info(f"Processing file: {args.input_faa}")

    # Load tokenizer
    tokenizer = AutoTokenizer.from_pretrained(esmt12_dir)

    # Step 1: Load and tokenize data
    step_start_time = time.time()
    logging.info("Step 1: Loading and tokenizing data...")

    # Call data loading and padding function
    data, seq_store_pos = load_and_tokenize_data(
        input_faa=args.input_faa,
        tokenizer=tokenizer,
        sample_source_path=args.negative_sample_path,
        max_length=args.sequence_length,
        batch_size=args.batch_sizes
    )

    # Extract encoded sequences and padded indices
    encoded_sequences = data["encoded_sequences"]
    padded_indices = data["padded_indices"]

    # Log data stats
    num_original_sequences = len(seq_store_pos)
    num_padded_sequences = len(padded_indices)
    total_sequences = len(encoded_sequences)

    step1_time = time.time() - step_start_time
    logging.info("Data loading and tokenization complete.")
    logging.info(f"Original sequences: {num_original_sequences}, Padded sequences: {num_padded_sequences}, Total sequences: {total_sequences}")
    logging.info(f"Step 1 completed in {step1_time:.2f} seconds.")
    
    # Step 2: Load model and extract embeddings
    step_start_time = time.time()

    # Skip step if output already exists
    if os.path.exists(results_file) and os.path.exists(positive_fasta_file) and os.path.exists(negative_fasta_file):
        logging.info("All output files already exist. Skipping model loading and feature extraction.")
        if os.path.exists(embeddings_tensor_path):
            logging.info(f"Loading embeddings tensor from {embeddings_tensor_path}")
            embeddings_tensor = torch.load(embeddings_tensor_path)
        else:
            raise FileNotFoundError(f"Embeddings tensor file not found at {embeddings_tensor_path}.")
    else:
        # Step 2: Load model and extract features
        logging.info("Step 2: Loading model and extracting features...")
        pretrain_path = os.path.join(esmt12_dir, 'model.safetensors')
        model = LRVMForClf(esmt12_dir, pretrain_path, 2)

        embeddings_tensor = extract_features(
            model=model,
            encoded_sequences=encoded_sequences,
            step_size=args.batch_sizes,
            device=args.device,
            padded_indices=padded_indices,
        )

        # Save extracted feature embeddings
        torch.save(embeddings_tensor, embeddings_tensor_path)
        logging.info(f"Feature embeddings saved to {embeddings_tensor_path}")

    step2_time = time.time() - step_start_time
    logging.info(f"Step 2 completed in {step2_time:.2f} seconds.")
    
    # Step 3: Classification
    step_start_time = time.time()
    logging.info("Step 3: Classifying sequences...")

    model_mlp = SimpleClassifier(input_dim=480, output_dim=2, hidden_dim=960, nhead=8, num_layers=1, dropout=0.4)
    load_model(model_mlp, args.weights, strict=True)

    classification_results = classify_genes(
        model_mlp=model_mlp,
        embeddings_tensor=embeddings_tensor,
        input_fasta_path=args.input_faa,
        device=args.device,
        threshold=args.threshold,
        padded_indices=padded_indices
    )

    positive_ids = classification_results["positive_ids"]
    negative_ids = classification_results["negative_ids"]

    logging.info(f"Classification completed. Total positive IDs: {len(positive_ids)}, Total negative IDs: {len(negative_ids)}")
    step3_time = time.time() - step_start_time
    logging.info(f"Step 3 completed in {step3_time:.2f} seconds.")

    # Final metrics logging
    elapsed_time = time.time() - overall_start_time
    final_memory = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024

    if torch.cuda.is_available():
        gpu_memory = torch.cuda.max_memory_allocated() / 1024 / 1024
    else:
        gpu_memory = 0

    metrics = {
        "input_file": args.input_faa,
        "num_sequences": len(seq_store_pos),
        "step1_time": round(step1_time, 2),
        "step2_time": round(step2_time, 2),
        "step3_time": round(step3_time, 2),
        "total_time": round(elapsed_time, 2),
        "initial_memory_MB": round(initial_memory, 2),
        "final_memory_MB": round(final_memory, 2),
        "gpu_memory_MB": round(gpu_memory, 2)
    }

    metrics_path = os.path.join(tmp_dir, f"{file_name}_runtime_metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    logging.info(f"🟢 Runtime metrics saved to {metrics_path}")

    return

if __name__ == "__main__":
    main_light()