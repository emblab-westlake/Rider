import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib_venn import venn3


RIDER_FILE = "Rider_results_gt200.txt"
LUCAPROT_FILE = "Luca_filtered_label1_len_gt200.csv"
PALMSCAN_FILE = "palmscan_results.tsv"

OUTPUT_DIR = Path("final_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FILTERED_PALMSCAN_OUT = OUTPUT_DIR / "filtered_palmscan_proteins_lg20_len200.txt"
VENN_OUT = OUTPUT_DIR / "venn_diagram_filtered_proteins_rider_lg200_high_quality.pdf"
CONSENSUS_OUT = OUTPUT_DIR / "consensus_protein_ids_rider_ge2.txt"


def clean_protein_id(protein_id):
    protein_id = protein_id.strip()
    protein_id = protein_id.lstrip(">")
    protein_id = re.sub(r"\s+", "", protein_id)
    return protein_id


def load_rider_ids(file_path):
    ids = set()
    with open(file_path, "r") as f:
        for line in f:
            if line.strip():
                raw_id = line.strip().split("\t")[0]
                ids.add(clean_protein_id(raw_id))
    return ids


def load_lucaprot_ids(file_path):
    ids = set()
    with open(file_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "protein_id" in row:
                ids.add(clean_protein_id(row["protein_id"]))
    return ids


def load_filtered_palmscan_ids(file_path, out_path):
    ids = set()
    with open(file_path, "r") as f_in, open(out_path, "w") as f_out:
        header = f_in.readline()
        for line in f_in:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            try:
                score = float(parts[1])
                aa_length = int(parts[-1])
            except ValueError:
                continue

            if score > 20 and aa_length >= 200:
                pid = clean_protein_id(parts[0])
                ids.add(pid)
                f_out.write(pid + "\n")

    return ids


def save_id_list(id_set, out_path):
    with open(out_path, "w") as f:
        for pid in sorted(id_set):
            f.write(pid + "\n")


def main():
    rider_ids = load_rider_ids(RIDER_FILE)
    lucaprot_ids = load_lucaprot_ids(LUCAPROT_FILE)
    palmscan_ids = load_filtered_palmscan_ids(PALMSCAN_FILE, FILTERED_PALMSCAN_OUT)

    print(f"PalmScan (filtered proteins): {len(palmscan_ids)}")
    print(f"Rider (proteins): {len(rider_ids)}")
    print(f"LucaProt (proteins): {len(lucaprot_ids)}")

    print(f"PalmScan ∩ Rider: {len(palmscan_ids & rider_ids)}")
    print(f"PalmScan ∩ LucaProt: {len(palmscan_ids & lucaprot_ids)}")
    print(f"Rider ∩ LucaProt: {len(rider_ids & lucaprot_ids)}")
    print(f"All three intersect: {len(rider_ids & lucaprot_ids & palmscan_ids)}")

    plt.figure(figsize=(8, 8))
    venn3(
        [rider_ids, lucaprot_ids, palmscan_ids],
        ("Rider", "LucaProt", "PalmScan (Score>20 & Len≥200)"),
    )
    plt.title("Venn Diagram: Rider vs LucaProt vs PalmScan (Protein-Level)", fontsize=18)
    plt.tight_layout()
    plt.savefig(VENN_OUT, format="pdf")
    plt.close()

    print(f"Venn diagram saved to: {VENN_OUT}")

    rider_supported_ids = rider_ids & (lucaprot_ids | palmscan_ids)
    print(f"Protein IDs identified by Rider AND at least one other tool: {len(rider_supported_ids)}")

    save_id_list(rider_supported_ids, CONSENSUS_OUT)
    print(f"ID list saved to: {CONSENSUS_OUT}")


if __name__ == "__main__":
    main()