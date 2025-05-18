import os
import argparse

def extract_complete_sequences(input_faa, output_faa):
    """
    从 .faa 文件中提取以 'M' 开头并以 '*' 结尾的完整蛋白质序列，并去掉结尾的 '*'。

    :param input_faa: 输入的 .faa 文件路径
    :param output_faa: 输出的 .faa 文件路径，用于存储完整的蛋白质序列（结尾不含 '*'）
    """
    with open(input_faa, 'r') as infile, open(output_faa, 'w') as outfile:
        current_sequence = []
        sequence_id = None
        for line in infile:
            if line.startswith('>'):  # 序列ID行
                if current_sequence:
                    # 检查之前的序列是否完整
                    protein_sequence = ''.join(current_sequence).strip()
                    if protein_sequence.startswith('M') and protein_sequence.endswith('*'):
                        # 去掉结尾的 '*' 并写入输出文件
                        outfile.write(f"{sequence_id}{protein_sequence[:-1]}\n")
                # 更新序列ID并准备接收新序列
                sequence_id = line
                current_sequence = []
            else:
                # 收集序列行，去除多余的空白符
                current_sequence.append(line.strip())

        # 检查最后一个序列
        if current_sequence:
            protein_sequence = ''.join(current_sequence).strip()
            if protein_sequence.startswith('M') and protein_sequence.endswith('*'):
                outfile.write(f"{sequence_id}{protein_sequence[:-1]}\n")


def process_directory(input_dir, output_dir):
    """
    遍历指定目录下的所有 .faa 文件，提取完整的蛋白质序列并保存到输出目录。

    :param input_dir: 包含 .faa 文件的输入目录
    :param output_dir: 保存提取结果的输出目录
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 遍历输入目录中的所有 .faa 文件
    for filename in os.listdir(input_dir):
        if filename.endswith(".faa"):
            input_faa_file = os.path.join(input_dir, filename)
            output_faa_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_complete_protein.faa")

            print(f"正在处理文件: {input_faa_file}")
            extract_complete_sequences(input_faa_file, output_faa_file)
            print(f"提取结果已保存到: {output_faa_file}")


def parse_arguments():
    """
    解析命令行参数并返回相关选项。

    :return: 命令行参数对象
    """
    parser = argparse.ArgumentParser(description="提取 .faa 文件中的完整蛋白质序列")
    
    parser.add_argument("-a", "--input_faa",
                        type=str,
                        help="输入的 .faa 文件路径")
    
    parser.add_argument("-i", "--input_dir",
                        type=str,
                        help="包含 .faa 文件的输入目录路径")
    
    parser.add_argument("-u", "--output_faa",
                        type=str,
                        help="输出的 .faa 文件路径，用于存储提取的完整序列")
    
    parser.add_argument("-o", "--output_dir",
                        type=str,
                        default='./',
                        help="输出目录路径 (默认当前目录)")
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # 如果提供了输入目录，则处理目录中的 .faa 文件
    if args.input_dir:
        process_directory(args.input_dir, args.output_dir)
    # 否则，处理单个 .faa 文件
    elif args.input_faa and args.output_faa:
        extract_complete_sequences(args.input_faa, args.output_faa)
    else:
        print("请提供输入文件或输入目录。使用 -h 查看帮助信息。")