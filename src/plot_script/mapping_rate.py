#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
from pathlib import Path

# =============================================================================
# 1. NC-style global settings
# =============================================================================
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["font.size"] = 10

color_scatter = "#4D4D4D"
color_line = "#D62728"
color_box_edge = "#333333"
color_strip = "#1F77B4"

# =============================================================================
# 2. Input / Output
# =============================================================================
sample_suffix = "luca50"
# sample_suffix = "AS"
# sample_suffix = "Soil"
# sample_suffix = "INF"
# sample_suffix = "IBD9"

input_file = Path(
    f"final_mapping_cov_results_{sample_suffix}.txt"
)

output_pdf = Path(
    f"Analysis_rdssdb_{sample_suffix}.pdf"
)

print(f"Processing: {sample_suffix}")

# =============================================================================
# 3. Data loading
# =============================================================================
try:
    data = pd.read_csv(input_file, sep="\t", header=None)

    if data.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns, got {data.shape[1]}")

    data = data.iloc[:, :3]
    data.columns = ["Full_ID", "Score", "Mapping_Rate"]

    data["Score"] = pd.to_numeric(data["Score"], errors="coerce")
    data["Mapping_Rate"] = pd.to_numeric(data["Mapping_Rate"], errors="coerce")
    data = data.dropna(subset=["Score", "Mapping_Rate"]).copy()

except FileNotFoundError:
    print("⚠️ File not found. Using mock data for demonstration.")
    np.random.seed(42)
    n = 300
    s = np.random.uniform(60, 140, n)
    r = (s * 0.6) + np.random.normal(0, 8, n)
    data = pd.DataFrame({
        "Score": s,
        "Mapping_Rate": np.clip(r, 0, 100)
    })

if data.empty:
    raise RuntimeError("No valid data available for plotting.")

# =============================================================================
# 4. Binning
# =============================================================================
bin_step = 10
min_score = int(np.floor(data["Score"].min() / bin_step)) * bin_step
max_score = int(np.ceil(data["Score"].max() / bin_step)) * bin_step
bins = np.arange(min_score, max_score + bin_step, bin_step)
labels = [f"{i}-{i + bin_step}" for i in bins[:-1]]

data["Score_Interval"] = pd.cut(
    data["Score"],
    bins=bins,
    labels=labels,
    right=False,
    include_lowest=True
)

data = data.dropna(subset=["Score_Interval"]).copy()

# =============================================================================
# 5. Plotting
# =============================================================================
sns.set_style("white")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# --- Panel A: regression plot ---
sns.regplot(
    data=data,
    x="Score",
    y="Mapping_Rate",
    ax=ax1,
    scatter_kws={
        "s": 20,
        "alpha": 0.5,
        "color": color_scatter,
        "edgecolor": "none",
    },
    line_kws={
        "linewidth": 2,
        "color": color_line,
    },
    ci=95,
)

slope, intercept, r_value, p_value, std_err = stats.linregress(
    data["Score"], data["Mapping_Rate"]
)
r_sq = r_value ** 2

if p_value < 0.001:
    p_text = r"$P < 0.001$"
else:
    p_text = f"$P = {p_value:.3f}$"

stats_text = f"Pearson $r = {r_value:.2f}$\n$R^2 = {r_sq:.2f}$\n{p_text}"

ax1.text(
    0.05,
    0.95,
    stats_text,
    transform=ax1.transAxes,
    va="top",
    fontsize=9,
    linespacing=1.5,
)

ax1.set_xlabel("Score", fontsize=11, fontweight="bold")
ax1.set_ylabel("Mapping Rate (%)", fontsize=11, fontweight="bold")
ax1.set_title("a", loc="left", fontsize=14, fontweight="bold", pad=10)

# --- Panel B: boxplot + stripplot ---
sns.boxplot(
    data=data,
    x="Score_Interval",
    y="Mapping_Rate",
    ax=ax2,
    color="white",
    linewidth=1.2,
    width=0.6,
    showfliers=False,
    boxprops=dict(edgecolor=color_box_edge),
    whiskerprops=dict(color=color_box_edge),
    capprops=dict(color=color_box_edge),
    medianprops=dict(color=color_line, linewidth=1.5),
)

sns.stripplot(
    data=data,
    x="Score_Interval",
    y="Mapping_Rate",
    ax=ax2,
    color=color_strip,
    size=3,
    alpha=0.6,
    jitter=0.2,
    edgecolor="none",
)

ax2.set_xlabel("Score Interval", fontsize=11, fontweight="bold")
ax2.set_ylabel("Mapping Rate (%)", fontsize=11, fontweight="bold")
ax2.set_title("b", loc="left", fontsize=14, fontweight="bold", pad=10)
ax2.tick_params(axis="x", rotation=45)

# =============================================================================
# 6. Style cleanup
# =============================================================================
for ax in (ax1, ax2):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 5))
    ax.spines["bottom"].set_position(("outward", 5))

plt.tight_layout()
plt.subplots_adjust(wspace=0.25)

# =============================================================================
# 7. Save
# =============================================================================
plt.savefig(output_pdf, format="pdf", bbox_inches="tight", transparent=True)
plt.close()

print(f"✅ Figure saved to: {output_pdf}")