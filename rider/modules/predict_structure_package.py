"""
------------------------------------------
Script Name: predict_structure_package.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - predict rdrp protein structure
Notes:
    - Enjoy yourself!
------------------------------------------
"""
import os
import time
import torch
import concurrent.futures
from Bio import SeqIO
from transformers import EsmForProteinFolding
from tqdm import tqdm

# Function to load the model
def load_structure_model(model_path="/home/gaoyang/Rider/esmflod_v1", gpu_id=2):
    """
    Load the ESMFold model for protein structure prediction.

    :param model_path: Path to the model weights
    :param gpu_id: GPU ID to use
    :return: Loaded model
    """
    try:
        # os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        print(f"Using GPU: {os.environ['CUDA_VISIBLE_DEVICES']}")
        
        # Load the model
        model = EsmForProteinFolding.from_pretrained(model_path)
        config = model.config
        print(f"Initial chunk size: {config.esmfold_config.trunk.chunk_size}")
        config.esmfold_config.trunk.chunk_size = 64
        model = EsmForProteinFolding.from_pretrained(model_path, config=config)
        print(f"Modified chunk size: {model.config.esmfold_config.trunk.chunk_size}")
        model = model.eval().cuda()  # Move model to GPU
        return model
    except Exception as e:
        print(f"Error loading model: {e}")
        raise


# Function to generate a PDB file
def generate_pdb(record, model, output_dir, sequence_length):
    """
    Generate a PDB file from a protein sequence.

    :param record: A Bio.SeqIO record object containing sequence information
    :param model: Loaded ESMFold model
    :param output_dir: Output directory for the PDB file
    :param sequence_length: Length to truncate the sequence
    :return: Path to the generated PDB file and elapsed time
    """
    sequence_id = record.id
    sequence = str(record.seq)[:sequence_length]
    # print(f"Processing {sequence_id}, truncated sequence length: {len(sequence)}")

    # Check for invalid character 'X' or 'x' in the sequence
    if 'x' in sequence.lower():
        # print(f"Error processing {sequence_id}: Invalid character in the sequence: 'x'")
        return None, 0

    start_time = time.time()

    try:
        with torch.no_grad():
            output = model.infer_pdb(sequence)
    except Exception as e:
        print(f"Error processing sequence {sequence_id}: {e}")
        return None, 0

    os.makedirs(output_dir, exist_ok=True)
    pdb_file = os.path.join(output_dir, f"{sequence_id}.pdb")

    with open(pdb_file, 'w') as f:
        f.write(output)

    elapsed_time = time.time() - start_time
    # print(f"Generated {pdb_file}, time taken: {elapsed_time:.2f} seconds")

    return pdb_file, elapsed_time

# Batch processing function
def process_faa_files(input_dir, sequence_length, model, max_workers=1, suffix="_Rider_predicted_RNA_Virus_potential_candidates.faa"):
    """
    Batch process .faa files in the intermediate_dir of each subdirectory to generate corresponding PDB files.

    :param input_dir: Path to the root input directory
    :param sequence_length: Length to truncate each sequence
    :param model: Loaded ESMFold model
    :param max_workers: Maximum number of threads to use
    :param suffix: Suffix used to match .faa files
    """
    # Iterate over subdirectories in the root input directory
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)

        # Ensure it is a directory
        if os.path.isdir(subdir_path):
            # intermediate_dir is set to the subdirectory path
            intermediate_dir = subdir_path

            # Check if intermediate_dir exists
            if not os.path.exists(intermediate_dir):
                print(f"Intermediate directory not found: {intermediate_dir}")
                continue

            # Find all .faa files in intermediate_dir that match the given suffix
            faa_files = [
                os.path.join(intermediate_dir, file)
                for file in os.listdir(intermediate_dir)
                if file.endswith(suffix)
            ]

            # If no matching .faa files are found, skip this subdirectory
            if not faa_files:
                print(f"No matching .faa files found in {intermediate_dir}.")
                continue

            for input_faa_file in faa_files:
                print(f"Processing file: {input_faa_file}")

                # Construct output directory path
                base_name = os.path.basename(input_faa_file).replace(suffix, "")
                output_dir = os.path.join(os.path.dirname(input_faa_file), f"candidate_{sequence_length}_pdb")
                os.makedirs(output_dir, exist_ok=True)

                # Parse sequences from the .faa file
                records = list(SeqIO.parse(input_faa_file, "fasta"))

                # Filter out records that already have PDB files
                records_to_process = [record for record in records if not os.path.exists(os.path.join(output_dir, f"{record.id}.pdb"))]

                # Use tqdm to show progress bar with time estimate
                with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = [executor.submit(generate_pdb, record, model, output_dir, sequence_length) for record in records_to_process]

                    # Wrap thread pool with tqdm for progress tracking
                    for future in tqdm(concurrent.futures.as_completed(futures), total=len(records_to_process), desc="Processing sequences", unit="seq", colour="blue"):
                        pdb_file, elapsed_time = future.result()
                        if pdb_file:
                            print(f"Task {pdb_file} completed, time: {elapsed_time:.2f} sec")
