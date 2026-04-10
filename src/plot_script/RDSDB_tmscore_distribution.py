#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# 1. Illustrator Compatibility Settings
# =============================================================================
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "sans-serif"

# If needed, you can prefer a specific font:
# plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]

# =============================================================================
# 2. File Paths
# =============================================================================
FILE_PATH_POSITIVE = Path(
    "/root/gaoyang/rider_experiment/structure_benchmark/rider_RDSDB_distribution/rdrp_structure_db_analysis/rider_rdrp_stats_tmscore_top10.tsv"
)
FILE_PATH_NEGATIVE = Path(
    "/root/gaoyang/rider_experiment/structure_benchmark/rider_RDSDB_distribution/rdrp_structure_db_analysis/uniref90_negative_stats_with_tmscore_top10.tsv"
)

OUTPUT_DIR = Path(
    "/root/gaoyang/rider_experiment/structure_benchmark/rider_RDSDB_distribution/rdrp_structure_db_analysis/"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PDF = OUTPUT_DIR / "Rider_structure_database_tmscore_distribution_curve_compare.pdf"
OUTPUT_SVG = OUTPUT_DIR / "Rider_structure_database_tmscore_distribution_curve_compare.svg"

# =============================================================================
# 3. Data Processing Function
# =============================================================================
def load_and_process(file_path, group_name):
    """
    Read the TSV, split the TopK_tTMscores column by ';', explode into rows,
    and return only TMscore and Group columns.
    """
    try:
        df = pd.read_csv(file_path, sep="\t", on_bad_lines="skip", engine="python")

        if "TopK_tTMscores" not in df.columns:
            raise KeyError(f"'TopK_tTMscores' column not found in {file_path}")

        df["TopK_tTMscores"] = df["TopK_tTMscores"].astype(str).str.split(";")
        df_exploded = df.explode("TopK_tTMscores", ignore_index=True)

        df_exploded["TMscore"] = pd.to_numeric(
            df_exploded["TopK_tTMscores"].astype(str).str.strip(),
            errors="coerce"
        )

        df_exploded = df_exploded.dropna(subset=["TMscore"]).copy()
        df_exploded["Group"] = group_name

        return df_exploded[["TMscore", "Group"]]

    except Exception as e:
        print(f"Error processing {group_name} ({file_path}): {e}")
        return pd.DataFrame(columns=["TMscore", "Group"])


# =============================================================================
# 4. Load Data
# =============================================================================
print("Loading and processing data...")

df_pos = load_and_process(FILE_PATH_POSITIVE, "Positive")
df_neg = load_and_process(FILE_PATH_NEGATIVE, "Negative")

df_final = pd.concat([df_pos, df_neg], ignore_index=True)

print(f"Data loaded. Total rows for plotting: {len(df_final)}")
print(df_final.head())

if df_final.empty:
    raise RuntimeError("No valid data found for plotting.")


# =============================================================================
# 5. Visualization
# =============================================================================
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))

sns.kdeplot(
    data=df_final,
    x="TMscore",
    hue="Group",
    fill=True,
    alpha=0.5,
    linewidth=2.5,
    palette={"Positive": "purple", "Negative": "gold"},
    common_norm=False,
)

plt.title("TMscore Distribution Curve", fontsize=20, fontweight="bold", pad=20)
plt.xlabel("TM Score", fontsize=18)
plt.ylabel("Density", fontsize=18)
plt.tick_params(axis="both", which="major", labelsize=14)

sns.despine()
plt.tight_layout()

# =============================================================================
# 6. Save Output
# =============================================================================
plt.savefig(OUTPUT_PDF, format="pdf", transparent=True, bbox_inches="tight")
plt.savefig(OUTPUT_SVG, format="svg", transparent=True, bbox_inches="tight")
plt.close()

print("Plot saved successfully to:")
print(f"1. {OUTPUT_PDF}")
print(f"2. {OUTPUT_SVG}")