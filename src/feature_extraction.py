"""
------------------------------------------
Script Name: feature_extraction.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - Extract feature from aa sequence using protein language model
Notes:
    - Enjoy yourself!
------------------------------------------
"""

from tqdm import tqdm
import torch
import numpy as np


def extract_features(
    model,
    encoded_sequences,
    step_size=256,
    device="cuda",
    padded_indices=None
):
    """
    Extract feature embeddings from sequences using a Transformer model,
    and optionally filter out padded samples based on their indices.

    Args:
        model (torch.nn.Module): The deep learning model used for feature extraction.
        encoded_sequences (Dataset): Encoded sequences containing input tensors (input_ids).
        step_size (int): Batch size for inference.
        device (str): Device to use (default is "cuda").
        padded_indices (list, optional): Indices of padded samples to be excluded.

    Returns:
        torch.Tensor: Feature embedding tensor for valid samples (excluding padded ones).
    """
    print(f"Start sequence features extracting...")

    # Initialize variables
    embeddings = []
    max_size = len(encoded_sequences["input_ids"])
    start = 0

    # Set device
    device = torch.device(device)
    model = model.to(device)
    model.eval()

    # Extract features per batch
    with torch.no_grad():
        with tqdm(total=max_size, desc="Extracting features", unit="seq", colour="green") as pbar:
            while start < max_size:
                # Build current batch
                if start + step_size <= max_size:
                    input_id = encoded_sequences["input_ids"][start:start + step_size].to(device)
                else:
                    input_id = encoded_sequences["input_ids"][start:max_size].to(device)

                # Forward pass through the model
                outputs = model(input_id)
                features = outputs.features.cpu().numpy()  # Extract embeddings
                embeddings.append(features)

                # Free memory
                del input_id
                torch.cuda.empty_cache()

                # Update progress
                processed = min(step_size, max_size - start)
                start += processed
                pbar.update(processed)

    # Stack all embedding features
    embeddings_stack = np.vstack(embeddings)

    # Convert to tensor
    embeddings_tensor = torch.tensor(embeddings_stack)

    print(f"Feature extraction completed. Total sequences processed: {len(embeddings_tensor)}")
    return embeddings_tensor