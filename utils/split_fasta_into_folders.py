import os
import shutil
from pathlib import Path
import argparse

def distribute_files(source_dir, target_dirs, num_groups):
    all_files = sorted(Path(source_dir).glob("group_*.fasta"))
    num_files = len(all_files)
    files_per_group = num_files // num_groups

    for i, target_dir in enumerate(target_dirs):
        os.makedirs(target_dir, exist_ok=True)
        start_index = i * files_per_group
        end_index = start_index + files_per_group if i < num_groups - 1 else num_files
        for file in all_files[start_index:end_index]:
            shutil.move(str(file), target_dir)

def main():
    parser = argparse.ArgumentParser(description="Distribute FASTA files into multiple directories.")
    parser.add_argument("source_dir", help="Path to the source directory containing FASTA files.")
    parser.add_argument("target_dirs", nargs='+', help="List of target directories to distribute files into.")
    parser.add_argument("--num_groups", type=int, default=3, help="Number of groups to divide the files into.")

    args = parser.parse_args()
    distribute_files(args.source_dir, args.target_dirs, args.num_groups)

if __name__ == "__main__":
    main()

#python split_fasta_into_folders.py 
# path/to/source_directory 
# path/to/target_directory_1 
# path/to/target_directory_2 
# path/to/target_directory_3 
# --num_groups 3