"""
------------------------------------------
Script Name: mmcluster.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - mmseq cluster modual
Dependencies: 
    - Python 3.9 or higher
    - pandas >= 1.3.0
    - numpy >= 1.21.0
Notes:
    - Enjoy yourself!
------------------------------------------
"""
import os
import subprocess
import logging


def mmseqs2_clustering(
    known_rdrp_fasta_path, 
    predicted_rdrp_fasta_path, 
    tmp_dir, 
    file_name
):
    """
    Perform MMseqs2 clustering on known and predicted RDRP sequences.

    Args:
        known_rdrp_fasta_path (str): Path to the FASTA file containing known RDRP sequences.
        predicted_rdrp_fasta_path (str): Path to the FASTA file containing predicted RNA virus sequences.
        tmp_dir (str): Path to the temporary directory.
        file_name (str): Base name for output files.

    Returns:
        cluster_output_tsv (str): Path to the clustering result TSV file.
    """
    logging.info("Starting MMseqs2 clustering...")

    # Create temporary directories
    tmp_mmseq_dir = os.path.join(tmp_dir, "mmseqs")
    tmp_mmseq_tmp_dir = os.path.join(tmp_mmseq_dir, "tmp")
    os.makedirs(tmp_mmseq_dir, exist_ok=True)
    os.makedirs(tmp_mmseq_tmp_dir, exist_ok=True)

    # Define merged FASTA output path
    merged_rdrp_fasta_path = os.path.join(tmp_mmseq_dir, f"{file_name}_merged_rdrp.fasta")

    # Step 1: Merge FASTA files
    with open(merged_rdrp_fasta_path, "w") as outfile:
        # Write known RDRP sequences
        with open(known_rdrp_fasta_path, "r") as infile:
            outfile.write(infile.read())
        # Write predicted sequences
        with open(predicted_rdrp_fasta_path, "r") as infile:
            outfile.write(infile.read())
    logging.info(f"Merged fasta file created: {merged_rdrp_fasta_path}")

    # Step 2: Create MMseqs2 database
    mmseqs_db_path = os.path.join(tmp_mmseq_dir, file_name + ".mmseq2.db")
    mmseqs_createdb_cmd = [
        "mmseqs", "createdb",
        merged_rdrp_fasta_path,  # Input: merged RDRP FASTA
        mmseqs_db_path           # Output: MMseqs2 database
    ]
    subprocess.call(mmseqs_createdb_cmd)
    logging.info(f"MMseqs2 database created: {mmseqs_db_path}")

    # Step 3: Run MMseqs2 clustering
    mmseqs_cluster_db_path = os.path.join(tmp_mmseq_dir, file_name + ".mmseq2.cluster.db")
    mmseqs_cluster_cmd = [
        "mmseqs", "cluster",
        mmseqs_db_path,           # Input database
        mmseqs_cluster_db_path,   # Output clustered database
        tmp_mmseq_tmp_dir,        # Temporary directory
        "--min-seq-id", "0.9",    # Minimum sequence identity
        "-c", "0.333",            # Coverage threshold
        "-e", "0.1",              # E-value threshold
        "--cov-mode", "1",        # Coverage mode
        "--cluster-mode", "2"     # Clustering mode
    ]
    subprocess.call(mmseqs_cluster_cmd)
    logging.info(f"MMseqs2 clustering completed: {mmseqs_cluster_db_path}")

    # Step 4: Generate clustering output in TSV format
    cluster_output_tsv = os.path.join(tmp_dir, file_name + "_cluster_results.tsv")
    mmseqs_createtsv_cmd = [
        "mmseqs", "createtsv",
        mmseqs_db_path,           # Input database
        mmseqs_db_path,           # Target database
        mmseqs_cluster_db_path,   # Clustered database
        cluster_output_tsv        # Output TSV file
    ]
    subprocess.run(mmseqs_createtsv_cmd, check=True)
    logging.info(f"MMseqs2 clustering results saved to: {cluster_output_tsv}")

    # Step 5 (optional): Parse the clustering result to identify known genes
    # known_genes = set()
    # with open(cluster_output_tsv, "r") as f:
    #     for line in f:
    #         query, target, *_ = line.strip().split("\t")
    #         # If query and target are different, query is considered known
    #         if query != target:
    #             known_genes.add(query)
    # logging.info(f"Number of known genes identified: {len(known_genes)}")

    return cluster_output_tsv