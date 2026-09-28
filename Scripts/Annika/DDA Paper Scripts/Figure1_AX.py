import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


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
input_file = os.path.join(
    topdir,
    "Combined_Fragment_Assignment_Colored.xlsx"
)
output_folder = os.path.join(topdir, "AX")

datasets = ["AX"]

dataset_display_names = {
    "AX": "IQ-X AcquireX",
}

plot_datasets = [
    dataset_display_names.get(name, name)
    for name in datasets
]

dataset_colors = {
    "AX": "#ca3c97",
}

bar_colors = [
    dataset_colors[name]
    for name in datasets
]


# Calculate the confirmed percentage from the same confirmed and rejected
# ID populations used in Figure 1.
confirmed_df = pd.read_excel(
    input_file,
    sheet_name="Sheet1"
)

rejected_df = pd.read_excel(
    input_file,
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


true_values = [
    True,
    "TRUE",
    "True",
    "true",
    1,
    "1",
]

confirmed_df = confirmed_df[
    confirmed_df["Dataset"].isin(datasets)
    & confirmed_df["Overall Match Status"].isin(true_values)
].drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Lipid Match Key",
    ]
)

rejected_df = rejected_df[
    rejected_df["Dataset"].isin(datasets)
    & ~rejected_df["Overall Match Status"].isin(true_values)
].drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Lipid Match Key",
    ]
)


confirmed_counts = (
    confirmed_df
    .groupby("Dataset")
    .size()
    .reindex(datasets, fill_value=0)
)

rejected_counts = (
    rejected_df
    .groupby("Dataset")
    .size()
    .reindex(datasets, fill_value=0)
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


# Match Figure 1's overall dimensions/style.
fig = plt.figure(
    figsize=(6.5, 4.875),
    layout="constrained",
)

grid = fig.add_gridspec(
    nrows=1,
    ncols=3,
)

confirmed_ax = fig.add_subplot(grid[0, 0])
rejected_ax = fig.add_subplot(
    grid[0, 1],
    sharey=confirmed_ax,
)
percent_ax = fig.add_subplot(grid[0, 2])

axes = [
    confirmed_ax,
    rejected_ax,
    percent_ax,
]

x_positions = np.arange(
    len(plot_datasets)
)


# A: confirmed ID counts
confirmed_bars = confirmed_ax.bar(
    x_positions,
    confirmed_counts.values,
    color=bar_colors,
    edgecolor="black",
)

confirmed_ax.bar_label(
    confirmed_bars,
    labels=[
        str(value)
        for value in confirmed_counts.values
    ],
    padding=3,
    fontsize=8,
)

confirmed_ax.set_ylabel("Lipid IDs")
confirmed_ax.set_title(
    "Confirmed Lipid IDs",
    y=1.02
)


# B: rejected ID counts
rejected_bars = rejected_ax.bar(
    x_positions,
    rejected_counts.values,
    color=bar_colors,
    edgecolor="black",
)

rejected_ax.bar_label(
    rejected_bars,
    labels=[
        str(value)
        for value in rejected_counts.values
    ],
    padding=3,
    fontsize=8,
)

rejected_ax.set_ylabel("Lipid IDs")
rejected_ax.set_title(
    "Rejected Lipid IDs",
    y=1.02
)


# Match the same y-axis scaling logic used in Figure 1.
shared_count_limit = max(
    confirmed_counts.max(),
    rejected_counts.max(),
) * 1.12

confirmed_ax.set_ylim(
    0,
    shared_count_limit,
)


# C: percentage confirmed
percent_bars = percent_ax.bar(
    x_positions,
    true_positive_summary[
        "Percent Confirmed"
    ],
    color=bar_colors,
    edgecolor="black",
)

percent_ax.bar_label(
    percent_bars,
    labels=[
        f"{value:.1f}%"
        for value in true_positive_summary[
            "Percent Confirmed"
        ]
    ],
    padding=3,
    fontsize=8,
)

percent_ax.set_ylabel(
    "Confirmed IDs (%)"
)

percent_ax.set_title(
    "Percentage of IDs\nConfirmed",
    y=0.98
)

percent_ax.set_ylim(
    0,
    100,
)

percent_ax.yaxis.set_major_formatter(
    PercentFormatter(
        xmax=100,
        decimals=0,
    )
)


# Match Figure 1 axis formatting.
for ax in axes:
    ax.set_xticks(
        x_positions,
        labels=plot_datasets,
        rotation=0,
        ha="center",
    )

    ax.margins(x=0.35)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)


# Match Figure 1 panel-label placement.
for panel_label, ax in zip(
    ["A", "B", "C"],
    axes,
):
    ax.text(
        0.02,
        0.98,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="top",
    )


fig.savefig(
    os.path.join(
        output_folder,
        "Figure 1 AcquireX.png",
    ),
    dpi=600,
)

fig.savefig(
    os.path.join(
        output_folder,
        "Figure 1 AcquireX.pdf",
    ),
)


print(
    "\nConfirmed/rejected percentage summary:"
)

print(
    true_positive_summary.to_string()
)


plt.show(block=True)