"""
------------------------------------------
Script Name: predict_pipline.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
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
import torch
from torch import nn
from transformers import AutoTokenizer, EsmModel
from transformers.utils import ModelOutput
from Bio import SeqIO
from safetensors.torch import load_model
from datasets import Dataset
import argparse
import logging
from rider.utils.download import ensure_data_exists
from rider.modules.predict_structure_package import load_structure_model, process_faa_files
from rider.modules.mmcluster import mmseqs2_clustering  
from rider.modules.classification import classify_genes
from rider.modules.feature_extraction import extract_features
from rider.modules.structure_align import Structure_aligned
from rider.modules.filter_prob_taxonomy import process_fixed_prob_out_folders
from rider.modules.parse_result import process_final_result  

# logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
# get workdir
script_dir = os.path.dirname(os.path.abspath(__file__))
# get root dir
project_root = os.path.dirname(script_dir) 
# get negative calibration
default_negative_sample_path = os.path.join(script_dir, "databases", "false256.faa")
# Set default submodule & database paths
USER_DATA_DIR = os.path.expanduser("~/.rider_data")
# Set default path for submodule 
default_submodule_path = os.path.join(USER_DATA_DIR, "submodule")
# Set default path for the Foldseek RDRP structure database
default_rdrp_db_path = os.path.join(USER_DATA_DIR, "Rider_pdb_database", "database")
# Set default path for weight
default_weight_path=os.path.join(project_root, "checkpoint", "checkpoint-44000","model.safetensors") #Rider/checkpoint/checkpoint-102000/model.safetensors

# parse arg
def parse_args():
    parser = argparse.ArgumentParser(description="Protein sequence classification and structure prediction pipeline.")
    parser.add_argument("-i", "--input_faa", required=True, type=str, 
                        help="Path to the input FASTA file.")
    parser.add_argument("-w", "--weights", required=True, type=str, default=default_weight_path,
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
    parser.add_argument("--threshold", type=float, default=0.8, 
                        help="Threshold value for classification (default: 0.8).")
    parser.add_argument("--device", type=str, default="cuda", 
                        help="Device to use for computation (default: 'cuda').")
    parser.add_argument("--negative_sample_path", type=str, default=default_negative_sample_path, 
                        help="Path to the negative samples file.")
    parser.add_argument("--structure_align_enabled", action="store_true", 
                        help="Enable Foldseek alignment task.")
    parser.add_argument("--submodule_dir", type=str, default=default_submodule_path,
                        help="Path to the submodule directory (e.g., containing Foldseek, ESM2 weights, etc.)")
    parser.add_argument("--rdrp_structure_database", type=str, required=False, default=default_rdrp_db_path,
                        help="Path to Foldseek database.")
    parser.add_argument("--alignment-type", type=int, default=1, 
                        help="Foldseek alignment type parameter.")
    parser.add_argument("--prob_threshold", type=int, default=50, 
                        help="Homologous probability threshold default 50.")
    parser.add_argument("--top_n_mean_prob", type=int, default=1, 
                        help="Usually the at least top 1 query above the homo prob threshold, you can choose 2 or more to obtain more strict results.")
    
    args = parser.parse_args()

    # Validate that the database exists if structure alignment is enabled
    if args.structure_align_enabled:
        if not os.path.exists(args.rdrp_structure_database):
            parser.error(f"Foldseek database not found at {args.rdrp_structure_database}. "
                         f"Please download it and place it in the expected location or specify with --rdrp_structure_database.")

    
    if not os.path.exists(args.submodule_dir):
        print(f"⚠️ Submodule not found at {args.submodule_dir}")
        print("⏳ Attempting to download submodule...")
        ensure_data_exists(download_submodule=True)

    if args.structure_align_enabled and not os.path.exists(args.rdrp_structure_database):
        print(f"⚠️ Foldseek database not found at {args.rdrp_structure_database}")
        print("⏳ Attempting to download database...")
        ensure_data_exists(download_db=True)
    return args

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
    def __init__(self,
                 input_dim=480,
                 output_dim=2,
                 hidden_dim=960,
                 nhead=8,
                 num_layers=2,
                 dropout=0.3,
                 class_weights=None):
        """
        Simple classifier built on top of TransformerEncoder.
        - Adds optional class_weights buffer to be compatible with checkpoints that saved it.
        - Tries to construct TransformerEncoderLayer with batch_first=True if supported
          to avoid nested tensor warning and for better inference performance.
        """
        super(SimpleClassifier, self).__init__()

        # Try to create encoder layer with batch_first if supported (PyTorch >= 1.11)
        try:
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=input_dim,
                nhead=nhead,
                dim_feedforward=hidden_dim,
                dropout=dropout,
                batch_first=True  # better performance and avoids nested tensor warning
            )
            self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        except TypeError:
            # older PyTorch: batch_first not supported, fall back
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=input_dim,
                nhead=nhead,
                dim_feedforward=hidden_dim,
                dropout=dropout
            )
            self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        self.linear = nn.Linear(input_dim, output_dim)

        # Ensure a class_weights buffer exists so loading checkpoints that contain it won't fail.
        # If class_weights provided, register with those values; otherwise register a neutral placeholder.
        if class_weights is not None:
            # accept list/ndarray/tensor
            cw = torch.tensor(class_weights, dtype=torch.float32)
            self.register_buffer('class_weights', cw)
        else:
            # register placeholder for 2 classes by default (adjust if multi-class)
            # This buffer will be overwritten when loading a checkpoint that contains class_weights.
            self.register_buffer('class_weights', torch.tensor([1.0, 1.0], dtype=torch.float32))

    def forward(self, input_ids=None, labels=None):
        """
        input_ids: either [batch, seq_len, input_dim] or [batch, input_dim]
        returns ModelOutput with keys: logits, features, (loss if labels provided)
        """
        x = input_ids
        # If transformer_encoder was created with batch_first=True, it expects [batch, seq, dim]
        # If input is [batch, input_dim], we treat it as seq_len=1 embedding per sample (no attention), so handle both.
        if x.dim() == 2:
            # treat as already the token-level representation (no sequence dim)
            first_token_output = x
        elif x.dim() == 3:
            # standard: [batch, seq_len, dim]
            first_token_output = x[:, 0, :]
        else:
            raise ValueError(f"Unexpected transformer output/input dimensions: {x.dim()}")

        logits = self.linear(first_token_output)
        loss = None
        if labels is not None:
            # use stored class_weights if available (not requiring it to be None)
            if hasattr(self, 'class_weights') and self.class_weights is not None:
                weights = self.class_weights.to(logits.device)
                loss_fct = nn.CrossEntropyLoss(weight=weights)
            else:
                loss_fct = nn.CrossEntropyLoss()
            loss = loss_fct(logits, labels.long())

            return ModelOutput({'loss': loss, 'logits': logits, 'features': first_token_output})
        return ModelOutput({'logits': logits, 'features': first_token_output})

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

def main():
    # Get the path of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set the absolute path to foldseek_dir
    foldseek_dir = os.path.join(project_root, "submodule", "foldseek", "bin")
    
    # Dynamically update the PATH environment variable
    os.environ["PATH"] = f"{foldseek_dir}:{os.environ.get('PATH', '')}"

    # Optional: verify that PATH was updated successfully
    print(f"Updated PATH: {os.environ['PATH']}")
    
    args = parse_args()
    start_time = time.time()
    overall_start_time = start_time  # Record the start time of the entire script
    
    # Path configuration
    file_name = os.path.basename(args.input_faa)
    tmp_dir = os.path.join(args.output_dir, file_name, f"{file_name}_intermediate")
    tmp_tensor_dir = os.path.join(tmp_dir, "tmp_tensors")
    embeddings_tensor_path = os.path.join(tmp_tensor_dir, f"{file_name}_features_embeddings.pt")
    tmp_mmseq_dir = os.path.join(tmp_dir, "tmp_mmseqs")
    tmp_mmseq_tmp_dir = os.path.join(tmp_mmseq_dir, "tmp")
    predicted_rdrp_fasta_path = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_RNA_Virus_potential_candidates.faa")
    
    # Model paths
    esmfold_dir=os.path.join(project_root,"submodule", "esmfold_v1")
    esmt12_dir=os.path.join(project_root,"submodule", "esm2_t12_35M_UR50D")
    known_rdrp = os.path.join(script_dir, "databases", "NCBI_RNA_virus_refseq.fasta")

    # Output file paths
    results_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_results.txt")
    positive_fasta_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_RNA_Virus_potential_candidates.faa")
    negative_fasta_file = os.path.join(tmp_dir, f"{file_name}_Rider_predicted_nonRNA.faa")
    
    # Ensure necessary directories exist
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(tmp_tensor_dir, exist_ok=True)
    os.makedirs(tmp_mmseq_dir, exist_ok=True)
    os.makedirs(tmp_mmseq_tmp_dir, exist_ok=True)
    
    # Configure dynamic log file path
    log_file_path = os.path.join(tmp_dir, "rider_pipeline.log")
    os.makedirs(args.output_dir, exist_ok=True)  # Ensure output directory exists
    
    # Configure logging to both file and console
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove previous handlers to avoid duplication
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Log initialization
    logging.info(f"Logging to file: {log_file_path}")
    logging.info(f"Processing file: {args.input_faa}")

    # Load tokenizer
    tokenizer = AutoTokenizer.from_pretrained(esmt12_dir)

    # Step 1: Load and tokenize data
    step_start_time = time.time()
    logging.info("Step 1: Loading and tokenizing data...")

    if os.path.exists(embeddings_tensor_path):
        logging.info("Step 1 skipped: Tokenized data already exists.")
    else:
        data, seq_store_pos = load_and_tokenize_data(
            input_faa=args.input_faa,
            tokenizer=tokenizer,
            sample_source_path=args.negative_sample_path,
            max_length=args.sequence_length,
            batch_size=args.batch_sizes
        )

        encoded_sequences = data["encoded_sequences"]
        padded_indices = data["padded_indices"]

        num_original_sequences = len(seq_store_pos)
        num_padded_sequences = len(padded_indices)
        total_sequences = len(encoded_sequences)

        logging.info("Data loading and tokenization complete.")
        logging.info(f"Original sequences: {num_original_sequences}, Padded sequences: {num_padded_sequences}, Total sequences: {total_sequences}")
        logging.info(f"Step 1 completed in {time.time() - step_start_time:.2f} seconds.")

    # Step 2: Load model and extract feature embeddings
    step_start_time = time.time()
    if os.path.exists(results_file) and os.path.exists(positive_fasta_file) and os.path.exists(negative_fasta_file):
        logging.info("Step 2 skipped: All output files already exist.")
        if os.path.exists(embeddings_tensor_path):
            logging.info(f"Loading embeddings tensor from {embeddings_tensor_path}")
            embeddings_tensor = torch.load(embeddings_tensor_path)
        else:
            raise FileNotFoundError(f"Embeddings tensor file not found at {embeddings_tensor_path}.")
    else:
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

        torch.save(embeddings_tensor, embeddings_tensor_path)
        logging.info(f"Feature embeddings saved to {embeddings_tensor_path}")
    
    logging.info(f"Step 2 completed in {time.time() - step_start_time:.2f} seconds.")
    
    # Step 3: Classification
    step_start_time = time.time()
    if os.path.exists(results_file):
        logging.info("Step 3 skipped: Classification results already exist.")
    else:
        logging.info("Step 3: Classifying sequences...")
        model_mlp = SimpleClassifier(input_dim=480, output_dim=2, hidden_dim=960, nhead=8, num_layers=2, dropout=0.3)
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
        save_predicted_genes(
            output_dir=tmp_dir,
            file_name=file_name,
            seq_store_pos=seq_store_pos,
            target_positive_record_id=positive_ids,
            target_negative_record_id=negative_ids,
            input_faa=args.input_faa
        )
    
    logging.info(f"Step 3 completed in {time.time() - step_start_time:.2f} seconds.")

    # Step 4: Clustering to filter known RDRP
    step_start_time = time.time()
    cluster_output_tsv = os.path.join(tmp_dir, "rdrp_clustering.tsv")
    if os.path.exists(cluster_output_tsv):
        logging.info("Step 4 skipped: Clustering results already exist.")
    else:
        logging.info("Step 4: Clustering...")
        cluster_output_tsv = mmseqs2_clustering(
            known_rdrp_fasta_path=known_rdrp,
            predicted_rdrp_fasta_path=predicted_rdrp_fasta_path,
            tmp_dir=tmp_dir,
            file_name="rdrp_clustering"
        )
        logging.info(f"Clustering result saved at: {cluster_output_tsv}")
    logging.info(f"Step 4 completed in {time.time() - step_start_time:.2f} seconds.")

    # Step 5: Structure prediction (optional)
    step_start_time = time.time()
    if args.predict_structure:
        logging.info("Step 5: Predicting structures...")
        structure_model = load_structure_model(
            model_path=args.structure_model_path or esmfold_dir,
            gpu_id=1
        )
        process_faa_files(
            input_dir=os.path.join(args.output_dir, file_name),
            sequence_length=args.sequence_length,
            model=structure_model,
            max_workers=1
        )
        logging.info("Structure prediction completed.")
    else:
        logging.info("Structure prediction is disabled.")
    logging.info(f"Step 5 completed in {time.time() - step_start_time:.2f} seconds.")
    
    # Step 6: Structural alignment
    step_start_time = time.time()
    if args.structure_align_enabled:
        logging.info("Step 6: Running Foldseek for structural alignment...")
        foldseek_runner = Structure_aligned(
            input_dir=args.output_dir,
            database_dir=args.rdrp_structure_database,
            alignment_type=args.alignment_type,
            sequence_length=args.sequence_length,
        )
        foldseek_runner.foldseek_batch()
        logging.info("Foldseek alignment completed.")
    else:
        logging.info("Foldseek alignment is disabled.")
    logging.info(f"Step 6 completed in {time.time() - step_start_time:.2f} seconds.")

    # Step 7: Filter queries based on Foldseek results
    step_start_time = time.time()
    if args.predict_structure:
        logging.info("Step 7: Filtering queries based on Foldseek results...")
        process_fixed_prob_out_folders(
            input_root_dir=args.output_dir,
            alignment_type=args.alignment_type,
            n=args.top_n_mean_prob,
            prob_threshold=args.prob_threshold
        )
        logging.info("Filtering completed.")
    else:
        logging.info("Structure prediction is disabled; filtering queries is skipped.")
    logging.info(f"Step 7 completed in {time.time() - step_start_time:.2f} seconds.")

    # Step 8: Extract final candidate sequences
    step_start_time = time.time()
    if args.predict_structure:
        logging.info("Step 8: Extracting sequences for final candidates...")
        process_final_result(args.output_dir, sample_name=file_name)
        logging.info("Sequence extraction completed.")
    else:
        logging.info("Sequence extraction is disabled.")
    logging.info(f"Step 8 completed in {time.time() - step_start_time:.2f} seconds.")
        
    # Final execution time
    elapsed_time = time.time() - overall_start_time
    logging.info(f"Total execution time: {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
