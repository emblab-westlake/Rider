python Benchmark_construct_unparied_test_set.py \
  --pos_fasta low_similarity_benchmarking_positive_serratus140w_00_lt300.faa \
  --neg_fasta low_similarity_excluded_uniref50_202104_benchmarking_negative_00_filtered_cleaned_2143800seqs.faa \
  --output_prefix benchmarking_test_0similarity_serratus140w_lt300 \
  --pos_per_set 100 \
  --neg_per_set 4000 \
  --num_sets 3 \
  --seed 42