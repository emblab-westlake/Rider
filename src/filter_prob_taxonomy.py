"""
------------------------------------------
Script Name: filter_prob.py
Version: 1.0.0
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2024-11-23
Description: 
    - filter homo prob
Notes:
    - Enjoy yourself!
------------------------------------------
"""

import os
import pandas as pd

def process_m8_file(m8_file_path, output_dir, prob_threshold=57):
    """
    Process a .m8 file to find the highest probability (prob) and corresponding target
    for each query, and save the results.

    :param m8_file_path: Path to the .m8 file
    :param output_dir: Directory to save the output .txt file
    :param prob_threshold: Threshold to filter the probability (prob)
    """
    col_names = ['query', 'target', 'qca', 'tca', 'alntmscore', 'qtmscore', 
                 'ttmscore', 'u', 't', 'lddt', 'lddtfull', 'prob']
    
    # Read the .m8 file
    df = pd.read_csv(m8_file_path, sep='\t', names=col_names)

    # Initialize a list to store qualified (query, target, prob) tuples
    results = []

    # Group by 'query' and find the target with the highest prob for each query
    for query, group in df.groupby('query'):
        # Exclude rows where target is the same as query or target starts with 'Rider'
        filtered_group = group[(group['target'] != query) & (~group['target'].str.startswith('Rider'))]

        # Skip this query if no valid target remains
        if filtered_group.empty:
            continue

        # Select the row with the highest `prob`
        best_match = filtered_group.loc[filtered_group['prob'].idxmax()]

        # Retrieve the target and prob
        best_target = best_match['target']
        best_prob = best_match['prob']

        # Only keep queries with prob >= prob_threshold
        if best_prob >= prob_threshold:
            results.append((query, best_target, best_prob))

    # Print the number of queries that meet the criteria
    print(f"Number of queries that meet the criteria: {len(results)}")

    # If there are any valid queries, save them to a .txt file
    if results:
        output_file_path = os.path.join(output_dir, "final_result_taxonomy_id.txt")

        with open(output_file_path, 'w') as output_file:
            for query, target, prob in results:
                # Remove 'Rider_' prefix from query
                cleaned_query = query.replace("Rider_", "")
                output_file.write(f"{cleaned_query}\t{target}\t{prob:.2f}\n")

        print(f"Results saved to {output_file_path}")
    else:
        print("No queries met the criteria.")


def process_fixed_prob_out_folders(input_root_dir, alignment_type=1, n=1, prob_threshold=57):
    """
    Process the .m8 files under the _intermediate/prob_out folders in each subdirectory of input_root_dir.
    Filters queries based on the given probability threshold and saves the results as .txt files.

    :param input_root_dir: Root input directory containing subdirectories for each sample
    :param alignment_type: Alignment type parameter used to select the appropriate .m8 files
    :param n: Number of top prob values to average (currently unused in this function)
    :param prob_threshold: Threshold to filter probability (prob)
    """

    # Iterate over first-level subdirectories
    for subdir in os.listdir(input_root_dir):
        subdir_path = os.path.join(input_root_dir, subdir)

        # Ensure it's a directory
        if os.path.isdir(subdir_path):
            # Construct the _intermediate folder path
            intermidiate_dir = os.path.join(subdir_path, f"{subdir}_intermediate")
            # intermidiate_dir = os.path.join(subdir_path)

            # Ensure the _intermediate folder exists
            if os.path.exists(intermidiate_dir):
                # Construct the prob_out folder path
                prob_out_dir = os.path.join(intermidiate_dir, "prob_out")
                
                # Check if prob_out folder exists
                if os.path.exists(prob_out_dir):
                    # Look for .m8 files in the prob_out directory
                    for filename in os.listdir(prob_out_dir):
                        if filename.endswith(f"_aln_type{alignment_type}.m8"):
                            m8_file_path = os.path.join(prob_out_dir, filename)
                            print(f"Processing file: {m8_file_path}")
                            
                            # Save results in the same subdirectory
                            process_m8_file(m8_file_path, subdir_path, n=n, prob_threshold=prob_threshold)
                else:
                    print(f"Skipped: prob_out folder does not exist: {prob_out_dir}")
            else:
                print(f"Skipped: _intermediate folder does not exist: {intermidiate_dir}")


if __name__ == "__main__":
    # Example input directory
    input_dir = "/usr/commondata/public/gaoyang/software/rider/project/yuanlin_AS"
    
    # Call the batch processing function
    process_fixed_prob_out_folders(input_dir, alignment_type=1, n=1, prob_threshold=50)
