cd /root/gaoyang/Rider/baseline/scripts
for i in /root/gaoyang/Rider/baseline/test_set/serratus140w_100_4000/*.csv; do
  echo "Processing file: $i"
  CUDA_VISIBLE_DEVICES=1 python 01_generate_tokenID.py \
    --csv-path "$i" \
    --output-pt /root/gaoyang/Rider/baseline/test_set/serratus140w_100_4000/tokenID/$(basename "${i%.csv}_token_ids.pt") \
    --model-path /root/gaoyang/Rider/esm2_t12_35M_UR50D \
    --max-length 1024
done

#/root/gaoyang/Rider/baseline/test_set/benchmark_25_25_10000_hard