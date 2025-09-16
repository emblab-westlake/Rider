"""
------------------------------------------
Script Name: parse_final_result.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - parse final result to generate candidate aa sequence into a fasta file.
Notes:
    - Enjoy yourself!
------------------------------------------
"""
import os

def parse_faa_file(faa_file):
    """
    Parse a .faa file and return a dictionary where the key is the sequence ID 
    and the value is the sequence string.

    :param faa_file: Path to the .faa file
    :return: Dictionary of sequence ID -> sequence content
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(faa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # If a new sequence ID is encountered, save the previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Extract new sequence ID and remove 'Rider_' prefix
                raw_id = line[1:].split()[0]  # Take the part after '>' up to the first space
                current_id = raw_id.replace("Rider_", "")
                current_seq = []  # Reset the sequence buffer
            else:
                # Accumulate sequence lines
                current_seq.append(line)
        
        # Save the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def extract_sequences_from_faa(faa_file, final_result_file, output_file):
    """
    Extract sequences from a .faa file based on query IDs listed in final_result.txt, 
    and save the matching sequences and their probabilities into an output file.

    :param faa_file: Path to the .faa file
    :param final_result_file: Path to the final_result.txt file (contains query ID and prob)
    :param output_file: Path to save the extracted results
    """
    # Read final_result.txt to get query IDs and their probs
    query_prob_dict = {}
    with open(final_result_file, 'r') as f:
        for line in f:
            query, prob = line.strip().split()
            query_prob_dict[query] = prob

    # Parse the .faa file to get all sequences
    sequences = parse_faa_file(faa_file)

    # Open output file to write matched sequences
    found_any = False
    with open(output_file, 'w') as out_f:
        for query, prob in query_prob_dict.items():
            if query in sequences:
                # Write ID, sequence, and prob
                out_f.write(f"{query}\t{sequences[query]}\t{prob}\n")
                found_any = True
            else:
                print(f"Warning: ID {query} not found in {faa_file}")
        
        # Ensure the file is created even if no matches are found
        if not found_any:
            out_f.write("")  # Create an empty file


def process_final_result(output_dir, sample_name=None):
    """
    Traverse subdirectories under output_dir, process final_result.txt files,
    extract corresponding sequences from .faa files, and save them as final_candidate.txt.

    :param output_dir: Path to the output directory
    :param sample_name: Name of the specific sample to process (optional)
    """
    for subdir in os.listdir(output_dir):
        # Skip directories that do not match the sample_name (if provided)
        if sample_name and subdir != sample_name:
            continue

        subdir_path = os.path.join(output_dir, subdir)

        # Ensure it's a directory
        if os.path.isdir(subdir_path):
            # Construct the _intermediate folder path
            intermidiate_dir = os.path.join(subdir_path, f"{subdir}_intermediate")
            faa_file = os.path.join(intermidiate_dir, f"{subdir}_Rider_predicted_RNA_Virus_potential_candidates.faa")
            print(f"Parsing faa file path: {faa_file}")
            final_result_file = os.path.join(subdir_path, "final_result.txt")
            output_file = os.path.join(subdir_path, "final_candidate.txt")

            # Check existence of .faa file and final_result.txt
            if not os.path.exists(faa_file):
                print(f"Error: .faa file not found: {faa_file}")
                continue

            if not os.path.exists(final_result_file):
                print(f"Error: final_result.txt file not found: {final_result_file}")
                continue

            # Extract sequences and save
            extract_sequences_from_faa(faa_file, final_result_file, output_file)
            print(f"Final candidate sequences saved to {output_file}")