
cd /root/gaoyang/Rider/baseline/scripts
for i in /root/gaoyang/Rider/baseline/test_set/serratus140w_100_4000/*.csv; do
  echo "Processing file: $i"
  CUDA_VISIBLE_DEVICES=0 python 01_generate_embeddings.py \
    --input_csv $i \
    --output_pt /root/gaoyang/Rider/baseline/test_set/serratus140w_100_4000/ESM2_35M_embeddings/$(basename "${i%.csv}_embeddings.pt") \
    --esm_model_path /root/gaoyang/Rider/submodule/esm2_t12_35M_UR50D
done

#/root/gaoyang/Rider/baseline/test_set/benchmark_25_25_10000_hard