#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path


FONT_CONFIG = {
    "title_size": 18,
    "x_label_size": 16,
    "x_tick_size": 18,
    "y_label_size": 18,
    "y_tick_size": 16,
    "legend_size": 16,
    "legend_title_size": 12,
}

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

RIDER_FILE = "Rider_results_gt200.txt"
LUCAPROT_FILE = "Luca_filtered_label1_len_gt200.csv"
PALMSCAN_FILE = "palmscan_global_all_complete_proteins_with_length.tsv"

OUTPUT_DIR = Path("final_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FILTERED_PALMSCAN_FILE = OUTPUT_DIR / "filtered_palmscan_proteins_lg20_len200.txt"
OUTPUT_PDF = OUTPUT_DIR / "environment_comparison_gt200.pdf"


def count_environments(file_path):
    counts = {}
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            env = line.split("_")[0]
            counts[env] = counts.get(env, 0) + 1
    return counts


def filter_palmscan(file_path, out_path):
    filtered_ids = set()
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
                filtered_ids.add(parts[0])
                f_out.write(parts[0] + "\n")

    return filtered_ids


def main():
    filter_palmscan(PALMSCAN_FILE, FILTERED_PALMSCAN_FILE)

    env_counts1 = count_environments(RIDER_FILE)
    env_counts2 = count_environments(LUCAPROT_FILE)
    env_counts3 = count_environments(FILTERED_PALMSCAN_FILE)

    all_environments = sorted(
        set(env_counts1.keys()) | set(env_counts2.keys()) | set(env_counts3.keys())
    )

    df = pd.DataFrame({
        "Rider": [env_counts1.get(env, 0) for env in all_environments],
        "LucaProt": [env_counts2.get(env, 0) for env in all_environments],
        "Palmscan": [env_counts3.get(env, 0) for env in all_environments],
    }, index=all_environments)

    df["Total"] = df.sum(axis=1)
    df = df.sort_values(by="Total", ascending=False).drop(columns=["Total"])

    df_log = np.log10(df + 1)

    bar_width = 0.25
    x = np.arange(len(df_log))

    fig, ax = plt.subplots(figsize=(14, 7))

    ax.bar(x - bar_width, df_log["Rider"], width=bar_width, label="Rider", color="#4C72B0", zorder=3)
    ax.bar(x, df_log["LucaProt"], width=bar_width, label="LucaProt", color="#DD8452", zorder=3)
    ax.bar(x + bar_width, df_log["Palmscan"], width=bar_width, label="Palmscan", color="#55A868", zorder=3)

    ax.set_xlabel("Environment", fontsize=FONT_CONFIG["x_label_size"], labelpad=10)
    ax.set_ylabel("log10(Counts + 1)", fontsize=FONT_CONFIG["y_label_size"], labelpad=10)
    ax.set_title("Environment Counts Comparison", fontsize=FONT_CONFIG["title_size"], pad=20)

    ax.set_xticks(x)
    ax.set_xticklabels(
        df_log.index,
        rotation=45,
        ha="right",
        fontsize=FONT_CONFIG["x_tick_size"],
    )

    max_val = df_log.max().max()
    y_limit_int = int(np.ceil(max_val))
    y_ticks = range(0, y_limit_int + 1)
    y_tick_labels = [f"$10^{i}$" if i > 0 else "0" for i in y_ticks]

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=FONT_CONFIG["y_tick_size"])
    ax.set_ylim(0, y_limit_int + 0.5)
    ax.set_xlim(-0.5 - bar_width, len(df_log) - 1 + 0.5 + bar_width)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)

    ax.legend(
        title="Method",
        fontsize=FONT_CONFIG["legend_size"],
        title_fontsize=FONT_CONFIG["legend_title_size"],
        frameon=False,
        loc="upper right",
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_PDF, format="pdf", transparent=True, dpi=300)
    plt.close()

    print(f"Plot saved to: {OUTPUT_PDF}")
    print(f"Font config: {FONT_CONFIG}")


if __name__ == "__main__":
    main()