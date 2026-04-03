"""
------------------------------------------
Script Name: filter_prob_taxonomy.py
Version: 1.9.0 (Updated: Added TM-score filtering logic)
Author: Gaoyang Luo
Contact information: lgyjsnjhit@gmail.com
                     luogaoyang@westlake.edu.cn
                     gaoyang.luo@unsw.edu.au
Created Date: 2024-11-20
Last Modified Date: 2026-03-27
Description: 
    - Filter Foldseek results based on Bits score or TM-score.
    - Output 1: final_result_taxonomy_id.txt (ID format: SRR..._1;1-33)
    - Output 2: final_result_taxonomy_id_with_seq.txt (High Quality Filtered)
    - Output 3: final_result_taxonomy_id_with_seq_low_quality.txt (Filtered out results)
------------------------------------------
"""

import os
import pandas as pd

def load_raw_sequences(txt_path):
    """
    Reads a raw 3-column TXT file (No header).
    Format: ID  Sequence  Label
    Returns a dictionary {id: sequence}.
    """
    sequences = {}
    if not os.path.exists(txt_path):
        print(f"Warning: Sequence file not found: {txt_path}")
        return sequences

    with open(txt_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Split by whitespace (tab or space)
            parts = line.split()
            
            # Ensure there are at least 2 columns (ID and Sequence)
            if len(parts) >= 2:
                seq_id = parts[0]  # The ID in txt file (e.g., SRR5947895_k99_11266_2)
                seq_data = parts[1]
                sequences[seq_id] = seq_data
    
    return sequences

def process_m8_file(m8_file_path, output_dir, seq_file_path, prob_threshold=50, top_n=1, threshold_type=1):
    """
    Process a .m8 file. 
    threshold_type: 1=bits, 2=ttmscore, 3=qtmscore, 4=max(ttmscore, qtmscore)
    """
    col_names = [
        'query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'tstart', 'tend', 
        'evalue', 'lddt', 'ttmscore', 'bits', 'qtmscore'
    ]
    
    try:
        df = pd.read_csv(m8_file_path, sep='\t', names=col_names)
        
        numeric_cols = ['bits', 'qtmscore', 'ttmscore', 'lddt', 'evalue', 'qstart', 'qend']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            
        df = df.dropna(subset=['bits', 'qtmscore', 'ttmscore']) 
        
    except pd.errors.EmptyDataError:
        print(f"Skipped: Empty file {m8_file_path}")
        return

    # Calculate tmscore threshold
    tmscore_threshold = prob_threshold / 100.0

    # Add max_tmscore column for type 4
    if threshold_type == 4:
        df['max_tmscore'] = df[['ttmscore', 'qtmscore']].max(axis=1)

    results = []

    for query, group in df.groupby('query'):
        filtered_group = group[(group['target'] != query) & (~group['target'].str.startswith('Rider'))]

        if filtered_group.empty:
            continue

        # Sort based on threshold_type
        if threshold_type == 1:
            sorted_group = filtered_group.sort_values('bits', ascending=False)
        elif threshold_type == 2:
            sorted_group = filtered_group.sort_values('ttmscore', ascending=False)
        elif threshold_type == 3:
            sorted_group = filtered_group.sort_values('qtmscore', ascending=False)
        elif threshold_type == 4:
            sorted_group = filtered_group.sort_values('max_tmscore', ascending=False)
        else:
            sorted_group = filtered_group.sort_values('bits', ascending=False)

        top_n_hits = sorted_group.head(top_n)

        avg_bits = top_n_hits['bits'].mean()
        avg_qtm = top_n_hits['qtmscore'].mean()
        avg_ttm = top_n_hits['ttmscore'].mean()
        avg_lddt = top_n_hits['lddt'].mean()

        best_match = top_n_hits.iloc[0]
        best_target = best_match['target']
        best_evalue = best_match['evalue']
        best_qstart = int(best_match['qstart']) 
        best_qend = int(best_match['qend'])     

        # Filter based on selected threshold_type
        meets_criteria = False
        if threshold_type == 1 and avg_bits >= prob_threshold:
            meets_criteria = True
        elif threshold_type == 2 and avg_ttm >= tmscore_threshold:
            meets_criteria = True
        elif threshold_type == 3 and avg_qtm >= tmscore_threshold:
            meets_criteria = True
        elif threshold_type == 4 and max(avg_ttm, avg_qtm) >= tmscore_threshold:
            meets_criteria = True

        if meets_criteria:
            results.append((query, best_target, avg_bits, avg_qtm, avg_ttm, avg_lddt, best_evalue, best_qstart, best_qend))

    print(f"Number of queries that meet the criteria: {len(results)}")

    if results:
        output_file_path = os.path.join(output_dir, "final_result_taxonomy_id.txt")
        with open(output_file_path, 'w') as output_file:
            for item in results:
                query, target, bits, qtm, ttm, lddt, evalue = item[:7]
                cleaned_query_file1 = query.replace("Rider_", "")
                output_file.write(f"{cleaned_query_file1}\t{target}\t{bits:.0f}\t{qtm:.4f}\t{ttm:.4f}\t{lddt:.2f}\t{evalue:.2e}\n")
        print(f"Original results saved to {output_file_path}")

        seq_dict = load_raw_sequences(seq_file_path)
        
        output_target_path = os.path.join(output_dir, "final_result_taxonomy_id_with_seq.txt")
        output_rest_path = os.path.join(output_dir, "final_result_taxonomy_id_with_seq_low_quality.txt")
        
        count_target = 0
        count_rest = 0

        with open(output_target_path, 'w') as f_tgt, open(output_rest_path, 'w') as f_rst:
            for query, target, bits, qtm, ttm, lddt, evalue, qstart, qend in results:
                temp_query = query.replace("Rider_", "") 
                start_pos, end_pos = "NA", "NA"
                id_part_full = temp_query 
                
                if ";" in temp_query:
                    id_part_full, coords = temp_query.split(";", 1) 
                    if "-" in coords:
                        try:
                            start_pos, end_pos = coords.split("-", 1)
                        except ValueError:
                            pass 
                
                if "_" in id_part_full:
                    clean_id = id_part_full.rsplit("_", 1)[0]
                else:
                    clean_id = id_part_full

                sequence = seq_dict.get(clean_id, "Sequence_Not_Found")

                is_high_quality = False
                try:
                    s_pos = float(start_pos)
                    e_pos = float(end_pos)
                    length = abs(s_pos - e_pos)
                    
                    # Determine which score and thresholds to use for quality filtering
                    if threshold_type == 1:
                        score = bits
                        thresh1, thresh2 = 50, 60
                    elif threshold_type == 2:
                        score = ttm
                        thresh1, thresh2 = 0.5, 0.6
                    elif threshold_type == 3:
                        score = qtm
                        thresh1, thresh2 = 0.5, 0.6
                    elif threshold_type == 4:
                        score = max(ttm, qtm)
                        thresh1, thresh2 = 0.5, 0.6

                    if length >= 100 and score >= thresh1:
                        is_high_quality = True
                    elif 80 <= length <= 99 and score >= thresh2:
                        is_high_quality = True
                        
                except ValueError:
                    is_high_quality = False

                line_to_write = (
                    f"{clean_id}\t{target}\t{bits:.0f}\t{qtm:.4f}\t{ttm:.4f}\t{lddt:.2f}\t{evalue:.2e}\t"
                    f"{qstart}\t{qend}\t{start_pos}\t{end_pos}\t{sequence}\n"
                )

                if is_high_quality:
                    f_tgt.write(line_to_write)
                    count_target += 1
                else:
                    f_rst.write(line_to_write)
                    count_rest += 1
        
        print(f"Extended results processing complete:")
        print(f"  - High Quality ({count_target} lines): {output_target_path}")
        print(f"  - Low Quality  ({count_rest} lines): {output_rest_path}")

    else:
        print("No queries met the criteria.")


def process_fixed_prob_out_folders(input_root_dir, alignment_type=1, n=1, prob_threshold=50, threshold_type=1):
    """
    Process the .m8 files under the _intermediate/prob_out folders.
    """
    for subdir in os.listdir(input_root_dir):
        subdir_path = os.path.join(input_root_dir, subdir)

        if os.path.isdir(subdir_path):
            intermediate_dir = os.path.join(subdir_path, f"{subdir}_intermediate")

            if os.path.exists(intermediate_dir):
                seq_filename = f"{subdir}_Rider_predicted_results.txt" 
                seq_file_path = os.path.join(intermediate_dir, seq_filename)
                
                if not os.path.exists(seq_file_path):
                    print(f"Warning: Raw sequence file not found at expected path: {seq_file_path}")

                prob_out_dir = os.path.join(intermediate_dir, "prob_out")
                
                if os.path.exists(prob_out_dir):
                    for filename in os.listdir(prob_out_dir):
                        if filename.endswith(f"_aln_type{alignment_type}.m8"):
                            m8_file_path = os.path.join(prob_out_dir, filename)
                            print(f"Processing file: {m8_file_path} (n={n}, type={threshold_type})")
                            
                            process_m8_file(
                                m8_file_path, 
                                subdir_path, 
                                seq_file_path, 
                                prob_threshold=prob_threshold, 
                                top_n=n,
                                threshold_type=threshold_type
                            )


if __name__ == "__main__":
    input_dir = "/usr/commondata/public/gaoyang/20250117_store/global_case_adjustment2/g5/group_5.fasta"
    process_fixed_prob_out_folders(input_dir, alignment_type=1, n=1, prob_threshold=50, threshold_type=1)