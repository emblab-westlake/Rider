python Benchmark_construct_unparied_test_set.py \
  --pos_fasta low_similarity_benchmarking_positive_00_127343seqs.faa \
  --neg_fasta low_similarity_excluded_uniref50_202104_benchmarking_negative_00_filtered_cleaned_2143800seqs.faa \
  --other_fasta benchmark_3protein.fasta \
  --output_prefix benchmarking_test_0similarity_with_hard_neg3_proteins \
  --pos_per_set 25 \
  --neg_per_set 10000 \
  --num_sets 4 \
  --seed 42