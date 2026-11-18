"""Plot lipid class composition and unique lipid counts by acquisition method.

Classes use Figure 4's grouping of related Ontology subclasses. Accepted
IDs are selected using Overall Match Status. Following Figure 4, each distinct
grouped class / Metabolite name before '|' counts once per dataset, across duplicate
features, adducts, and polarities. Summed annotations are included: this analysis
requires no chain information. Counts describe identifications, not intensity.

Percentages are 100 * unique IDs in a class / all unique IDs in that dataset.
All classes are shown by default. --min-class-percent optionally selects classes
reaching that percentage in any dataset, without changing denominators. Missing
classes are zero; datasets with no eligible IDs display N/A.

Save SI_Class_Heatmap as PNG (600 dpi) and PDF, with panel A showing
raw counts and panel B showing percentages. Save only PNG/PDF figures.
Input can be Excel or CSV. Example:
python SI_Class_Heatmap.py --input Combined_Fragment_Assignment_Colored.xlsx --no-show
Use --all-ids explicitly for unvalidated annotations.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np
import pandas as pd

from Figure3 import simplify_lipid_class


INPUT_FILE = Path(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline\Combined_Fragment_Assignment_Colored.xlsx")
DATASETS = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
DATASET_LABELS = {"OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT",
                  "Stellar": "Stellar", "Stellarboth": "Stellar (\u00b1)"}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def prepare_lipids(source, datasets, all_ids=False):
    """Select IDs and deduplicate using Figure 4's grouped class/name identity."""
    required = {"Dataset", "Ontology", "Metabolite name"}
    if not all_ids:
        required.add("Overall Match Status")
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Missing columns: " + ", ".join(sorted(missing)))
    source = source.copy()
    for column in ["Dataset", "Ontology", "Metabolite name"]:
        source[column] = source[column].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False
        ).str.strip()
    source["Dataset"] = source["Dataset"].str.casefold().map(
        {dataset.casefold(): dataset for dataset in datasets}
    )
    eligible = source["Dataset"].notna() & source["Ontology"].ne("")
    if not all_ids:
        eligible &= source["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    source["Lipid Match Key"] = source["Metabolite name"].str.split("|", n=1).str[0].str.strip()
    eligible &= source["Lipid Match Key"].ne("")
    lipids = source.loc[eligible].copy()
    lipids["Lipid Class"] = lipids["Ontology"].apply(simplify_lipid_class)
    lipids = lipids.drop_duplicates(["Dataset", "Lipid Class", "Lipid Match Key"])
    if lipids.empty:
        raise ValueError("No eligible lipid IDs with class and name in the selected datasets.")
    return lipids


def class_tables(lipids, datasets, min_class_percent=0):
    """Return displayed counts/percentages and totals before row filtering."""
    if not np.isfinite(min_class_percent) or not 0 <= min_class_percent <= 100:
        raise ValueError("Minimum class percentage must be between 0 and 100.")
    counts = pd.crosstab(lipids["Lipid Class"], lipids["Dataset"]).reindex(
        columns=datasets, fill_value=0
    ).sort_index()
    denominators = counts.sum(axis=0).astype(int).rename("Total unique lipid IDs")
    percentages = counts.div(denominators.replace(0, np.nan), axis=1) * 100
    selected = percentages.ge(min_class_percent).any(axis=1)
    counts = counts.loc[selected].astype("Int64")
    counts.loc[:, denominators.eq(0)] = pd.NA
    return counts, denominators, percentages.loc[selected]


def plot_class_panel(ax, counts, percentages, denominators, population, raw):
    table = counts if raw else percentages
    values = table.astype(float).to_numpy()
    if not table.empty:
        cmap = plt.get_cmap("Blues").copy()
        cmap.set_bad("#E6E6E6")
        vmax = max(1, float(np.nanmax(values))) if raw else 100
        heatmap = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=0,
                            vmax=vmax, aspect="auto", interpolation="nearest")
        for row in range(len(table)):
            for column in range(len(table.columns)):
                value = values[row, column]
                label = "N/A"
                if np.isfinite(value):
                    label = f"{value:,.0f}" if raw else f"{value:.1f}%"
                ax.text(column, row, label, ha="center", va="center",
                        color="white" if value >= 0.6 * vmax else "black", fontsize=8)
        ax.set_xticks(np.arange(len(table.columns) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(table) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="minor", bottom=False, left=False)
        colorbar = ax.figure.colorbar(heatmap, ax=ax, pad=0.03)
        colorbar.set_label("Unique Lipid IDs" if raw else "Unique Lipid IDs (%)")
        if raw:
            colorbar.locator = MaxNLocator(integer=True)
            colorbar.update_ticks()
        else:
            colorbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    else:
        ax.set_xlim(-0.5, len(table.columns) - 0.5)
        ax.text(0.5, 0.5, "No classes meet the threshold", transform=ax.transAxes,
                ha="center", va="center")
    ax.set_xticks(np.arange(len(table.columns)), labels=[
        f"{DATASET_LABELS.get(dataset, dataset)}\nIDs = {denominators[dataset]:,}"
        for dataset in table.columns
    ], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(table)), labels=table.index)
    ax.set_xlabel("Acquisition Method", labelpad=12)
    ax.set_ylabel("Lipid Class")
    measure = "Counts" if raw else "Composition"
    ax.set_title(f"{population} Lipid {measure}", pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1)


def plot_class_heatmap(counts, percentages, denominators, population):
    """Combine raw counts (A) and percentages (B) in manuscript format."""
    fig, axes = plt.subplots(
        1, 2, figsize=(8, max(4.875, 0.22 * len(counts) + 1.8)),
        layout="constrained",
    )
    for ax, label, raw in zip(axes, "AB", [True, False]):
        plot_class_panel(ax, counts, percentages, denominators, population, raw)
        ax.text(0, 1.04, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", ha="left", va="bottom")
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--min-class-percent", type=float, default=0)
    parser.add_argument("--all-ids", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    if not np.isfinite(args.min_class_percent) or not 0 <= args.min_class_percent <= 100:
        parser.error("--min-class-percent must be between 0 and 100")
    datasets = list(dict.fromkeys(args.datasets))
    source = pd.read_csv(args.input) if args.input.suffix.lower() == ".csv" else pd.read_excel(args.input)
    lipids = prepare_lipids(source, datasets, args.all_ids)
    counts, denominators, percentages = class_tables(lipids, datasets, args.min_class_percent)
    output_dir = args.output_dir or args.input.parent / "Final Figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    population = "Annotated" if args.all_ids else "Accepted"
    fig = plot_class_heatmap(counts, percentages, denominators, population)
    for extension in ["png", "pdf"]:
        path = output_dir / f"SI_Class_Heatmap.{extension}"
        fig.savefig(path, dpi=600)
        print(f"Saved {path}")
    if args.no_show:
        plt.close(fig)
    print("\nUnique lipid IDs per dataset:")
    print(denominators.to_string())
    if not args.no_show:
        plt.show(block=True)


if __name__ == "__main__":
    main()
