#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


RIDER_FILE = "Rider_filtered_results_gt200.txt"
LUCAPROT_FILE = "Luca_filtered_label1_len_gt200.csv"
PALMSCAN_FILE = "palmscan_global_all_complete_proteins_with_length.tsv"

OUTPUT_DIR = Path("final_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PDF = OUTPUT_DIR / "protein_length_distribution_deparate_contigs_gt200_high_quality_updated_ge2.pdf"

MODE = "ge2"          # any, ge2, rider, luca, palm, custom
LENGTH_MAX = 4000
TRUNCATE_MODE = "filter"   # filter or truncate


def extract_contig_id(protein_id):
    if protein_id is None:
        return ""
    protein_id = str(protein_id).strip().lstrip(">")
    protein_id = re.sub(r"\s+", "", protein_id)
    return re.sub(r"_[0-9]+$", "", protein_id)


def load_rider_ids(file_path):
    ids = set()
    with open(file_path, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts:
                contig_id = extract_contig_id(parts[0])
                if contig_id:
                    ids.add(contig_id)
    return ids


def load_luca_ids(file_path):
    ids = set()
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            protein_id = row.get("protein_id") or row.get("proteinId") or row.get("id") or ""
            contig_id = extract_contig_id(protein_id)
            if contig_id:
                ids.add(contig_id)
    return ids


def load_palm_ids(file_path):
    ids = set()
    with open(file_path, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                score = float(parts[1])
                aa_len = int(parts[-1])
            except ValueError:
                continue

            if score > 20 and aa_len >= 200:
                contig_id = extract_contig_id(parts[0])
                if contig_id:
                    ids.add(contig_id)
    return ids


def get_consensus_ids(set1, set2, set3, mode):
    if mode == "ge2":
        return (set1 & set2) | (set1 & set3) | (set2 & set3)
    if mode == "rider":
        return set1
    if mode == "luca":
        return set2
    if mode == "palm":
        return set3
    if mode == "any":
        return set1 | set2 | set3
    if mode == "custom":
        return set1 | (set2 & set3)
    raise ValueError(f"Unknown MODE: {mode}")


def collect_rider_lengths(consensus_ids, file_path):
    lengths = []
    with open(file_path, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= 11:
                continue
            contig_id = extract_contig_id(parts[0])
            if contig_id in consensus_ids:
                aa_seq = parts[11]
                lengths.append(len(aa_seq))
    return lengths


def collect_luca_lengths(consensus_ids, file_path):
    lengths = []
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            protein_id = row.get("protein_id") or row.get("proteinId") or row.get("id") or ""
            contig_id = extract_contig_id(protein_id)
            if contig_id in consensus_ids:
                seq = row.get("seq") or row.get("sequence") or row.get("aa") or ""
                lengths.append(len(seq))
    return lengths


def collect_palm_lengths(consensus_ids, file_path):
    lengths = []
    with open(file_path, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                score = float(parts[1])
                aa_len = int(parts[-1])
            except ValueError:
                continue

            contig_id = extract_contig_id(parts[0])
            if contig_id in consensus_ids and score > 20 and aa_len >= 200:
                lengths.append(aa_len)
    return lengths


def apply_length_filter(lengths, max_len):
    return [x for x in lengths if x <= max_len]


def safe_median(arr):
    return int(np.median(arr)) if arr else None


def main():
    set1 = load_rider_ids(RIDER_FILE)
    set2 = load_luca_ids(LUCAPROT_FILE)
    set3 = load_palm_ids(PALMSCAN_FILE)

    consensus_ids = get_consensus_ids(set1, set2, set3, MODE)
    print(f"MODE={MODE}, consensus size: {len(consensus_ids)}")

    rider_lengths = collect_rider_lengths(consensus_ids, RIDER_FILE)
    luca_lengths = collect_luca_lengths(consensus_ids, LUCAPROT_FILE)
    palm_lengths = collect_palm_lengths(consensus_ids, PALMSCAN_FILE)

    print(f"Raw counts -> Rider: {len(rider_lengths)}, LucaProt: {len(luca_lengths)}, PalmScan: {len(palm_lengths)}")

    if TRUNCATE_MODE == "filter":
        rider_plot = apply_length_filter(rider_lengths, LENGTH_MAX)
        luca_plot = apply_length_filter(luca_lengths, LENGTH_MAX)
        palm_plot = apply_length_filter(palm_lengths, LENGTH_MAX)
    else:
        rider_plot = rider_lengths[:]
        luca_plot = luca_lengths[:]
        palm_plot = palm_lengths[:]

    print(
        f"After filter (<= {LENGTH_MAX}): "
        f"Rider {len(rider_plot)}, LucaProt {len(luca_plot)}, PalmScan {len(palm_plot)}"
    )

    plt.figure(figsize=(8, 6))
    colors = {
        "Rider": "#1f77b4",
        "LucaProt": "#ff7f0e",
        "PalmScan": "#2ca02c",
    }

    if rider_plot:
        sns.kdeplot(rider_plot, label="Rider", color=colors["Rider"], linewidth=2.0, bw_adjust=0.5)
    if luca_plot:
        sns.kdeplot(luca_plot, label="LucaProt", color=colors["LucaProt"], linewidth=2.0, bw_adjust=0.5)
    if palm_plot:
        sns.kdeplot(palm_plot, label="PalmScan", color=colors["PalmScan"], linewidth=2.0, bw_adjust=0.5)

    r_m = safe_median(rider_plot if TRUNCATE_MODE == "filter" else rider_lengths)
    l_m = safe_median(luca_plot if TRUNCATE_MODE == "filter" else luca_lengths)
    p_m = safe_median(palm_plot if TRUNCATE_MODE == "filter" else palm_lengths)

    if r_m is not None:
        plt.axvline(r_m, color=colors["Rider"], linestyle="--", alpha=0.6)
    if l_m is not None:
        plt.axvline(l_m, color=colors["LucaProt"], linestyle="--", alpha=0.6)
    if p_m is not None:
        plt.axvline(p_m, color=colors["PalmScan"], linestyle="--", alpha=0.6)

    if TRUNCATE_MODE == "truncate":
        plt.xlim(0, LENGTH_MAX)

    title_map = {
        "ge2": "Protein Length Distribution for Consensus Contigs (≥2 Tools)",
        "rider": "Protein Length Distribution for Rider-identified Contigs",
        "luca": "Protein Length Distribution for LucaProt-identified Contigs",
        "palm": "Protein Length Distribution for PalmScan-identified Contigs",
        "any": "Protein Length Distribution for Contigs Identified by Any Tool",
        "custom": "Protein Length Distribution (Custom Rule)",
    }

    plt.xlabel("Protein Length (aa)", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    plt.title(title_map.get(MODE, "Protein Length Distribution"), fontsize=16, fontweight="bold")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved to: {OUTPUT_PDF}")


if __name__ == "__main__":
    main()