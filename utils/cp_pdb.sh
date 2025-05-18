#!/bin/bash

# 输入文件路径
input_file="/home/gaoyang/Rider/plot_draw/IBD_RPKM/final_filtered_len_100_to_200.txt"

# 基础路径
base_dir="/usr/commondata/public/gaoyang/software/rider/RNA_virus_project/IBD_human_mt_sortmerna"
output_dir="$base_dir/all_result/pdb_gt200_all_mapping"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"

# 逐行读取输入文件
while read -r line; do
    # 获取第一列作为 ID（忽略后面的描述和数值）
    id=$(echo "$line" | awk '{print $1}')
    
    # 跳过空行
    [[ -z "$id" ]] && continue

    # 提取前缀，如 SRR5947972，从 ID 中分离出来
    prefix=$(echo "$id" | cut -d'_' -f1)

    # 构造完整路径
    pdb_path="${base_dir}/${prefix}.contigs.1K.fa_prodigalgv.faa/${prefix}.contigs.1K.fa_prodigalgv.faa_intermediate/candidate_1024_pdb/Rider_${id}.pdb"

    # 检查文件是否存在并复制
    if [ -f "$pdb_path" ]; then
        echo "✅ 正在复制: $pdb_path"
        cp "$pdb_path" "$output_dir/"
    else
        echo "⚠️ 未找到文件: $pdb_path"
    fi

done < "$input_file"

echo "🎉 所有复制任务完成！"