import argparse
from Bio import SeqIO

def filter_fasta(input_file, output_file, min_length=1000):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            if len(record.seq) > min_length:
                SeqIO.write(record, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences in a FASTA file by minimum length.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to the output FASTA file.")
    parser.add_argument("--min_length", type=int, default=1000, help="Minimum sequence length to retain (default: 1000).")

    args = parser.parse_args()
    filter_fasta(args.input_file, args.output_file, args.min_length)

if __name__ == "__main__":
    main()