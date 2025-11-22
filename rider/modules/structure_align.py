"""
------------------------------------------
Script Name: structure_align.py
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
import subprocess
import time
import logging

class Structure_aligned:
    def __init__(self, input_dir, database_dir, sequence_length, alignment_type=1):
        """
        Initialize the Structure_aligned class.

        :param input_dir: Root path to the input sample directories
        :param database_dir: Path to the Foldseek database
        :param alignment_type: Alignment type parameter (default is 1)
        :param sequence_length: Length used for truncating sequences
        """
        self.input_dir = input_dir
        self.database_dir = database_dir
        self.alignment_type = alignment_type
        self.sequence_length = sequence_length

    def run_foldseek(self, input_pdb_dir, output_file, tmp_dir):
        """
        Run Foldseek alignment on the specified PDB directory.

        :param input_pdb_dir: Directory containing sample PDB files
        :param output_file: Path to output alignment result file
        :param tmp_dir: Path for storing temporary files
        """
        foldseek_cmd = [
            "foldseek", "easy-search",
            input_pdb_dir,                # Sample PDB directory
            self.database_dir,            # Foldseek database directory
            output_file,                  # Output file path
            tmp_dir,                      # Temporary directory
            "--alignment-type", str(self.alignment_type)  # Alignment type parameter
        ]

        logging.info(f"Running Foldseek command: {' '.join(foldseek_cmd)}")

        # Record start time
        start_time = time.time()

        # Execute Foldseek command
        try:
            subprocess.run(foldseek_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Foldseek command failed: {e}")
            return None

        # Record elapsed time
        elapsed_time = time.time() - start_time
        logging.info(f"Foldseek completed in {elapsed_time:.2f} seconds")
        return output_file

    def foldseek_batch(self):
        """
        Run Foldseek alignment in batch mode for each sample directory.
        """
        for subdir in os.listdir(self.input_dir):
            subdir_path = os.path.join(self.input_dir, subdir)
        
            # Ensure it is a directory
            if os.path.isdir(subdir_path):
                # Construct the _intermediate directory path
                intermidiate_dir = os.path.join(subdir_path, f"{subdir}_intermediate")
                # print(f"intermidiate_dir: {intermidiate_dir}")
                # intermidiate_dir = subdir_path  # Alternative if structure differs

                # Check if _intermediate directory exists
                if os.path.exists(intermidiate_dir):
                    # Create prob_out directory
                    prob_out_dir = os.path.join(intermidiate_dir, "prob_out")
                    os.makedirs(prob_out_dir, exist_ok=True)

                    # Construct candidate_pdb directory path
                    candidate_pdb_dir = os.path.join(intermidiate_dir, f"candidate_{self.sequence_length}_pdb")
                    print(f"candidate_pdb_dir: {candidate_pdb_dir}")
                    if not os.path.exists(candidate_pdb_dir):
                        logging.warning(f"Skipping as candidate_pdb directory does not exist: {candidate_pdb_dir}")
                        continue
                        
                    # Construct output file path
                    sample_name = os.path.basename(subdir)
                    output_file = os.path.join(prob_out_dir, f"{sample_name}_aln_type{self.alignment_type}.m8")

                    # Check if output file already exists (Checkpoint)
                    if os.path.exists(output_file):
                        logging.info(f"Skipping Foldseek as output file already exists: {output_file}")
                        continue  # Skip current directory
                        
                    # Construct temporary directory path
                    tmp_dir = os.path.join(prob_out_dir, "tmp")
                    os.makedirs(tmp_dir, exist_ok=True)

                    # Run Foldseek alignment
                    self.run_foldseek(candidate_pdb_dir, output_file, tmp_dir)

                else:
                    logging.warning(f"Skipping as _intermediate directory does not exist: {intermidiate_dir}")