import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter, StrMethodFormatter


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
)

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

datasets = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
dataset_display_names = {
    "OTOT": "OT/OT",
    "OTIT": "OT/IT",
    "ITIT": "IT/IT",
    "Stellarboth": "Stellar (+/-)",
}
plot_datasets = [dataset_display_names.get(name, name) for name in datasets]

dataset_colors = {
    "OTOT": "#440154",
    "OTIT": "#414487",
    "ITIT": "#2A788E",
    "Stellar": "#22A884",
    "Stellarboth": "#7AD151",
}
bar_colors = [dataset_colors[name] for name in datasets]

# Calculate the confirmed percentage from the same confirmed and rejected
# ID populations used in Figure 1.
confirmed_df = pd.read_excel(
    "Combined_Fragment_Assignment_Colored.xlsx",
    sheet_name="Sheet1"
)
rejected_df = pd.read_excel(
    "Combined_Fragment_Assignment_Colored.xlsx",
    sheet_name="Sheet1"
)

for plot_df in [confirmed_df, rejected_df]:
    for column in ["Dataset", "Ontology"]:
        plot_df[column] = (
            plot_df[column]
            .astype(str)
            .str.replace("\u00a0", " ", regex=False)
            .str.strip()
        )
    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

true_values = [True, "TRUE", "True", "true", 1, "1"]
confirmed_df = confirmed_df[
    confirmed_df["Dataset"].isin(datasets)
    & confirmed_df["Overall Match Status"].isin(true_values)
].drop_duplicates(
    subset=["Dataset", "Ontology", "Lipid Match Key"]
)
rejected_df = rejected_df[
    rejected_df["Dataset"].isin(datasets)
    & ~rejected_df["Overall Match Status"].isin(true_values)
].drop_duplicates(
    subset=["Dataset", "Ontology", "Lipid Match Key"]
)

confirmed_counts = (
    confirmed_df.groupby("Dataset").size().reindex(datasets, fill_value=0)
)
rejected_counts = (
    rejected_df.groupby("Dataset").size().reindex(datasets, fill_value=0)
)

true_positive_summary = pd.DataFrame(
    {
        "Confirmed_IDs": confirmed_counts,
        "Rejected_IDs": rejected_counts,
    }
)
true_positive_summary["Total_IDs"] = (
    true_positive_summary["Confirmed_IDs"]
    + true_positive_summary["Rejected_IDs"]
)
true_positive_summary["Percent Confirmed"] = (
    true_positive_summary["Confirmed_IDs"]
    / true_positive_summary["Total_IDs"]
    * 100
)
true_positive_summary.index = plot_datasets

# Sum paired positive/negative scan counts within each replicate, then
# calculate the mean and sample standard deviation across the three runs.
polarity_scan_counts = {
    "OTOT": {
        "negative": [8801, 8793, 8839],
        "positive": [8665, 8673, 8659],
    },
    "OTIT": {
        "negative": [32098, 32283, 32331],
        "positive": [29386, 29226, 29304],
    },
    "ITIT": {
        "negative": [33792, 33739, 33392],
        "positive": [29958, 29955, 30040],
    },
    "Stellar": {
        "negative": [96948, 96671, 96152],
        "positive": [91040, 90562, 90072],
    }
}

scan_replicates = {
    dataset: (
        np.asarray(polarity_counts["negative"])
        + np.asarray(polarity_counts["positive"])
    )
    for dataset, polarity_counts in polarity_scan_counts.items()
}
scan_replicates["Stellarboth"] = np.asarray([71252, 70746, 70960])

scan_means = np.asarray(
    [scan_replicates[dataset].mean() for dataset in datasets]
)
scan_standard_deviations = np.asarray(
    [scan_replicates[dataset].std(ddof=1) for dataset in datasets]
)

fig = plt.figure(figsize=(6.5, 4.875), layout="constrained")
grid = fig.add_gridspec(nrows=2, ncols=2)

confirmed_ax = fig.add_subplot(grid[1, 0])
rejected_ax = fig.add_subplot(grid[0, 1], sharey=confirmed_ax)
scan_ax = fig.add_subplot(grid[0, 0])
percent_ax = fig.add_subplot(grid[1, 1])
axes = [scan_ax, rejected_ax, confirmed_ax, percent_ax]
x_positions = np.arange(len(plot_datasets))

# C: confirmed ID counts
confirmed_bars = confirmed_ax.bar(
    x_positions,
    confirmed_counts.values,
    color=bar_colors,
    edgecolor="black"
)
confirmed_ax.bar_label(
    confirmed_bars,
    labels=[str(value) for value in confirmed_counts.values],
    padding=3,
    fontsize=8
)
confirmed_ax.set_ylabel("Lipid IDs")
confirmed_ax.set_title("Confirmed Lipid IDs")

# B: rejected ID counts, sharing the same y-axis scale as panel C
rejected_bars = rejected_ax.bar(
    x_positions,
    rejected_counts.values,
    color=bar_colors,
    edgecolor="black"
)
rejected_ax.bar_label(
    rejected_bars,
    labels=[str(value) for value in rejected_counts.values],
    padding=3,
    fontsize=8
)
rejected_ax.set_ylabel("Lipid IDs")
rejected_ax.set_title("Rejected Lipid IDs")

shared_count_limit = max(confirmed_counts.max(), rejected_counts.max()) * 1.12
confirmed_ax.set_ylim(0, shared_count_limit)

# D: percentage confirmed
percent_bars = percent_ax.bar(
    x_positions,
    true_positive_summary["Percent Confirmed"],
    color=bar_colors,
    edgecolor="black"
)
percent_ax.bar_label(
    percent_bars,
    labels=[
        f"{value:.1f}%"
        for value in true_positive_summary["Percent Confirmed"]
    ],
    padding=3,
    fontsize=8
)
percent_ax.set_ylabel("Confirmed IDs (%)")
percent_ax.set_title("Percentage of IDs Confirmed")
percent_ax.set_ylim(0, 100)
percent_ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))

# A: mean scan count with sample-SD error bars
scan_ax.bar(
    x_positions,
    scan_means,
    yerr=scan_standard_deviations,
    capsize=5,
    color=bar_colors,
    edgecolor="black",
    error_kw={
        "elinewidth": 1.5,
        "capthick": 1.5,
        "ecolor": "black"
    }
)
for x_position, scan_mean, scan_standard_deviation in zip(
    x_positions,
    scan_means,
    scan_standard_deviations,
):
    scan_ax.annotate(
        f"{scan_mean:,.0f}",
        xy=(x_position, scan_mean + scan_standard_deviation),
        xytext=(0, 4),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=8,
    )
scan_ax.set_ylabel("Mean MS/MS Scans per Run")
scan_ax.set_title("Mean MS/MS Scan Count")
scan_ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
scan_ax.margins(y=0.14)

for ax in axes:
    ax.set_xticks(x_positions, labels=plot_datasets, rotation=45, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

for panel_label, ax in zip(["A", "B", "C", "D"], axes):
    ax.text(
        0.02,
        0.98,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="top"
    )

fig.savefig(
    "Figure 1.png",
    dpi=600
)
fig.savefig(
    "Figure 1.pdf"
)

print("\nConfirmed/rejected percentage summary:")
print(true_positive_summary.to_string())

plt.show(block=True)
