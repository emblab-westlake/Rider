"""
------------------------------------------
Script Name: filter_prob_taxonomy.py
Version: 1.1.0 (Updated for Top-N Averaging)
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                      luogaoyang@westlake.edu.cn
                      gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2026-02-12
Description: 
    - Filter Foldseek results based on Bits score.
    - Supports averaging scores of the Top N hits while retaining the Best Hit's target ID.
    - Output includes structural metrics (LDDT, TM-score) and E-value.
------------------------------------------
"""

import os
import pandas as pd

def process_m8_file(m8_file_path, output_dir, prob_threshold=57, top_n=1):
    """
    Process a .m8 file. For each query:
    1. Find valid targets (excluding self and 'Rider').
    2. Select the top N hits based on 'bits'.
    3. Calculate the AVERAGE scores (bits, tm, lddt) of these N hits.
    4. Keep the Target ID and E-value of the BEST hit (the #1 hit).
    5. Filter based on the Average Bits score.

    :param m8_file_path: Path to the .m8 file
    :param output_dir: Directory to save the output .txt file
    :param prob_threshold: Threshold to filter the probability (Average Bits)
    :param top_n: Number of top hits to average
    """
    col_names = [
        'query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'tstart', 'tend', 
        'evalue', 'lddt', 'ttmscore', 'bits', 'qtmscore'
    ]
    
    # Read the .m8 file
    try:
        df = pd.read_csv(m8_file_path, sep='\t', names=col_names)
    except pd.errors.EmptyDataError:
        print(f"Skipped: Empty file {m8_file_path}")
        return

    # Initialize a list to store results
    results = []

    # Group by 'query'
    # Note: Loop is used here for complex logic (Top N avg + Best Target ID) which is harder to vectorize cleanly
    for query, group in df.groupby('query'):
        # 1. Filter: Exclude self-matches and targets starting with 'Rider'
        filtered_group = group[(group['target'] != query) & (~group['target'].str.startswith('Rider'))]

        # Skip if no valid target remains
        if filtered_group.empty:
            continue

        # 2. Sort by bits descending to find the best hits
        sorted_group = filtered_group.sort_values('bits', ascending=False)

        # 3. Take Top N hits
        # If there are fewer than N hits, it will just take all of them.
        top_n_hits = sorted_group.head(top_n)

        # 4. Calculate Averages for scores
        avg_bits = top_n_hits['bits'].mean()
        avg_qtm = top_n_hits['qtmscore'].mean()
        avg_ttm = top_n_hits['ttmscore'].mean()
        avg_lddt = top_n_hits['lddt'].mean()

        # 5. Retrieve the Best Target info (from the #1 hit, which is the first row)
        best_match = top_n_hits.iloc[0]
        best_target = best_match['target']
        best_evalue = best_match['evalue'] # E-value usually shouldn't be averaged arithmetically

        # 6. Filter based on the AVERAGE bits score
        if avg_bits >= prob_threshold:
            results.append((query, best_target, avg_bits, avg_qtm, avg_ttm, avg_lddt, best_evalue))

    # Print the number of queries that meet the criteria
    print(f"Number of queries that meet the criteria (Top {top_n} Avg >= {prob_threshold}): {len(results)}")

    # If there are any valid queries, save them to a .txt file
    if results:
        output_file_path = os.path.join(output_dir, "final_result_taxonomy_id.txt")

        with open(output_file_path, 'w') as output_file:
            for query, target, bits, qtm, ttm, lddt, evalue  in results:
                # Remove 'Rider_' prefix from query
                cleaned_query = query.replace("Rider_", "")
                # Note: bits is formatted as .0f, but it's an average now, so it might have decimals.
                # If you want to keep it integer-like, use .0f, otherwise .2f
                output_file.write(f"{cleaned_query}\t{target}\t{bits:.0f}\t{qtm:.4f}\t{ttm:.4f}\t{lddt:.2f}\t{evalue:.2e}\n")

        print(f"Results saved to {output_file_path}")
    else:
        print("No queries met the criteria.")


def process_fixed_prob_out_folders(input_root_dir, alignment_type=1, n=1, prob_threshold=57):
    """
    Process the .m8 files under the _intermediate/prob_out folders.
    
    :param n: Number of top prob values to average. Passed to process_m8_file.
    """

    # Iterate over first-level subdirectories
    for subdir in os.listdir(input_root_dir):
        subdir_path = os.path.join(input_root_dir, subdir)

        # Ensure it's a directory
        if os.path.isdir(subdir_path):
            # Construct the _intermediate folder path
            intermediate_dir = os.path.join(subdir_path, f"{subdir}_intermediate")

            # Ensure the _intermediate folder exists
            if os.path.exists(intermediate_dir):
                # Construct the prob_out folder path
                prob_out_dir = os.path.join(intermediate_dir, "prob_out")
                
                # Check if prob_out folder exists
                if os.path.exists(prob_out_dir):
                    # Look for .m8 files in the prob_out directory
                    for filename in os.listdir(prob_out_dir):
                        if filename.endswith(f"_aln_type{alignment_type}.m8"):
                            m8_file_path = os.path.join(prob_out_dir, filename)
                            print(f"Processing file: {m8_file_path} (n={n})")
                            
                            # Save results in the same subdirectory
                            # Pass 'n' to the processing function
                            process_m8_file(m8_file_path, subdir_path, prob_threshold=prob_threshold, top_n=n)
                else:
                    # Silent skip or debug print
                    pass
            else:
                pass


if __name__ == "__main__":
    # Example input directory
    input_dir = "/usr/commondata/public/gaoyang/software/rider/project/yuanlin_AS"
    
    # Call the batch processing function
    # n=1: Same as original logic (Average of 1 is itself)
    # n=3: Average scores of top 3 hits, but keep Target ID of the #1 hit
    process_fixed_prob_out_folders(input_dir, alignment_type=1, n=1, prob_threshold=50)