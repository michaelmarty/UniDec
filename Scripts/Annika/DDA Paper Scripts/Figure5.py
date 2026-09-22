import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


# Make the repository importable when this script is run directly.
repository_root = Path(__file__).resolve().parents[3]
if str(repository_root) not in sys.path:
    sys.path.insert(0, str(repository_root))

from PublicScripts.Lipids.LipidFunctions import network_analysis  # noqa: E402


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

input_file = "Combined_Fragment_Assignment_Colored.xlsx"
output_basename = "Figure 5"

df = pd.read_excel(input_file, sheet_name="Sheet1")

required_columns = {
    "Dataset",
    "Overall Match Status",
    "Matched Required Fragments",
    "Metabolite name",
    "Ontology",
    "Polarity",
    "Average Mz",
    "MS/MS spectrum",
}
missing_columns = required_columns.difference(df.columns)
if missing_columns:
    raise ValueError(
        "The input workbook is missing required columns: "
        + ", ".join(sorted(missing_columns))
    )

def normalized_text(series):
    return (
        series
        .fillna("")
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
        .str.casefold()
    )


# Exclude Apex before calculating either density distribution or the networks.
dataset_text = normalized_text(df["Dataset"])
apex_count = int(dataset_text.eq("apex").sum())
df = df.loc[~dataset_text.eq("apex")].copy().reset_index(drop=True)
print(f"Excluded {apex_count} Apex IDs from analysis.")

true_values = [True, "TRUE", "True", "true", 1, "1"]
df["ID Status"] = np.where(
    df["Overall Match Status"].isin(true_values),
    "True",
    "False",
)

# Required-fragment matches are recorded as entries such as "FA1: Yes" and
# "FA2: No". Only an exact "Yes" is a successful match; results such as
# "Yes but Ref Failed" are failures. Ignore rejection flags and
# disqualifying-fragment results appended after the required-fragment results.
fragment_text = df["Matched Required Fragments"].fillna("").astype(str)
required_fragment_text = fragment_text.str.replace(
    r"\s+(?:FLAGGED:|DISQUALIFYING FRAGMENTS FOUND:).*$",
    "",
    regex=True,
)
df["Fragments Found"] = required_fragment_text.str.count(
    r":\s*Yes\s*(?=;|$)"
)
df["Fragments Missing"] = required_fragment_text.str.count(
    r":\s*(?:Yes but Ref Failed|No and Ref Failed|No|Not Found)\s*(?=;|$)"
)
df["Fragments Evaluated"] = df["Fragments Found"] + df["Fragments Missing"]
df["Fragments Matched (%)"] = np.where(
    df["Fragments Evaluated"] > 0,
    100 * df["Fragments Found"] / df["Fragments Evaluated"],
    np.nan,
)

# Calculate the same class-network strength used in Figure 4. Each ID is
# compared only with IDs from the same lipid class and polarity. Strength is
# the sum of its MS/MS spectral-similarity connection weights.
network_tolerance = 0.3
minimum_network_size = 3
df["_Figure Row"] = np.arange(len(df))
df["_Network Class"] = normalized_text(df["Ontology"])
df["_Network Polarity"] = normalized_text(df["Polarity"])
df["Class Strength"] = np.nan

analyzed_groups = 0
for _, network_input in df.groupby(
    ["_Network Class", "_Network Polarity"],
    sort=True,
):
    if len(network_input) < minimum_network_size:
        continue

    local_input = network_input.reset_index(drop=True)
    network_df, _ = network_analysis(
        local_input,
        mz_tol=network_tolerance,
        header="Class_",
        add_all_cols=True,
    )
    df.loc[
        network_df["_Figure Row"].astype(int).to_numpy(),
        "Class Strength",
    ] = network_df["Class_strength"].to_numpy()
    analyzed_groups += 1

df = df.drop(
    columns=["_Figure Row", "_Network Class", "_Network Polarity"]
)
df["Log10 Class Strength"] = np.log10(
    pd.to_numeric(df["Class Strength"], errors="coerce").clip(lower=0) + 1
)
print(
    f"Calculated class strength for {analyzed_groups} "
    f"class/polarity networks."
)

status_colors = {
    "True": "#D55E00",
    "False": "#0072B2",
}
status_labels = {
    "True": "Confirmed",
    "False": "Rejected",
}

fig, axes = plt.subplots(ncols=2, figsize=(6.5, 2.75))

strength_values = df["Log10 Class Strength"].dropna()
strength_upper = strength_values.max()
strength_bins = np.linspace(0, strength_upper, 25)

features = [
    (
        "Fragments Matched (%)",
        np.linspace(0, 100, 12),
        "Required Fragments Matched",
    ),
    (
        "Log10 Class Strength",
        strength_bins,
        "log10(1 + Class Strength)",
    ),
]

for ax, (column, bins, xlabel) in zip(axes, features):
    for status in ["False", "True"]:
        values = df.loc[df["ID Status"] == status, column].dropna()
        ax.hist(
            values,
            bins=bins,
            density=True,
            alpha=0.60,
            color=status_colors[status],
            edgecolor="black",
            linewidth=0.7,
            label=status_labels[status],
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_xlim(bins[0], bins[-1])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)

axes[0].set_title("Fragment Confirmation")
axes[1].set_title("Class-Network Strength")
axes[0].xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))

for panel_label, ax in zip(["A", "B"], axes):
    ax.text(
        0.02,
        1.04,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

fig.tight_layout()
fig.savefig(
    f"{output_basename}.png",
    dpi=600,
)
fig.savefig(
    f"{output_basename}.pdf",
)

print("\nIDs included in required-fragment panel:")
print(
    df.dropna(subset=["Fragments Matched (%)"])
    .groupby("ID Status")
    .size()
    .to_string()
)
print("\nIDs included in class-strength panel:")
print(
    df.dropna(subset=["Log10 Class Strength"])
    .groupby("ID Status")
    .size()
    .to_string()
)

plt.show(block=True)
