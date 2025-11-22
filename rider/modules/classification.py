"""
------------------------------------------
Script Name: classification.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - The main workflow of Rider
Usage: 
    - from src.classification import *
Dependencies: 
    - Python 3.9 or higher
    - pandas >= 1.3.0
    - numpy >= 1.21.0
Notes:
    - Enjoy yourself!
------------------------------------------
"""
import torch
import numpy as np
from Bio import SeqIO
import torch.nn.functional as F
from tqdm import tqdm

def classify_genes(
    model_mlp,
    embeddings_tensor,
    input_fasta_path,
    device,
    step_size=128,
    threshold=0.5,
    padded_indices=None,
):
    """
    Classify embeddings using a classification model while filtering out the indices of padded samples.

Args:
    model_mlp (torch.nn.Module): The MLP model used for classification.
    embeddings_tensor (torch.Tensor): Tensor of extracted positive sample embeddings.
    input_fasta_path (str): Path to the input FASTA file.
    device (torch.device): Device to run the model on (CPU or GPU).
    step_size (int): Batch size for inference.
    threshold (float): Threshold for classification.
    padded_indices (list, optional): Indices of padded samples to be filtered out.

Returns:
    dict: Classification results for positive and negative samples.
    """
    print(f"Start predicting RdRP genes...")

    # Set the model to evaluation mode
    model_mlp = model_mlp.to(device)
    model_mlp.eval()

    # Store classification results
    res_positive = []

    # Inference
    start = 0
    total_sequences = len(embeddings_tensor)
    
    with torch.no_grad():
        with tqdm(total=total_sequences, desc="Classifying sequences", unit="seq", colour="green") as pbar:
            while start < len(embeddings_tensor):
                if start + step_size <= len(embeddings_tensor):
                    input_id = embeddings_tensor[start:start + step_size].to(device)
                else:
                    input_id = embeddings_tensor[start:].to(device)

                # Model outputs logits
                outputs = model_mlp(input_id)
                logits = outputs.logits  # Assuming the model output contains logits
                logits_after_softmax = torch.softmax(logits, dim=1)
                pred = (logits_after_softmax[:, 1] > threshold).long().cpu().numpy()
                res_positive.extend(pred)

                # Free up memory
                del input_id
                torch.cuda.empty_cache()

                # start += step_size
                # Update progress bar
                processed = min(step_size, total_sequences - start)
                start += processed
                pbar.update(processed)  # Update progress bar

    print(f"Classification completed. Total sequences processed: {len(res_positive)}")

    # Parse FASTA sequences
    fasta_ids = []
    fasta_seqs = []
    for record in SeqIO.parse(input_fasta_path, 'fasta'):
        fasta_ids.append(str(record.id))
        fasta_seqs.append(str(record.seq))

    # Extract indices of positive and negative samples
    positive_indices = [i for i, label in enumerate(res_positive) if label == 1]
    negative_indices = [i for i, label in enumerate(res_positive) if label == 0]

    # If padded_indices are provided, filter out the padded samples
    if padded_indices is not None:
        positive_indices = [i for i in positive_indices if i not in padded_indices]
        negative_indices = [i for i in negative_indices if i not in padded_indices]

    # Extract IDs and sequences for positive and negative samples
    positive_ids = [fasta_ids[idx] for idx in positive_indices]
    positive_seqs = [fasta_seqs[idx] for idx in positive_indices]
    negative_ids = [fasta_ids[idx] for idx in negative_indices]
    negative_seqs = [fasta_seqs[idx] for idx in negative_indices]

    # Print results (only for original samples, excluding padded negatives)
    print(f"Total positive samples: {len(positive_ids)}")
    print(f"Total negative samples: {len(negative_ids)}")

    return {
        "positive_ids": positive_ids,
        "positive_seqs": positive_seqs,
        "negative_ids": negative_ids,
        "negative_seqs": negative_seqs,
    }