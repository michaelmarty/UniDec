"""Plot Figure 3: lipid-class coverage as a single panel."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


# Make the checkout's PublicScripts package available when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


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

input_file = "Combined_Fragment_Assignment_Colored.xlsx"
output_basename = "Figure3"
datasets = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]

# Muted sequential colors emphasize that coverage is ordered from one to five
# datasets. This is deliberately not viridis because viridis identifies
# individual datasets in the other figures.
coverage_colors = {
    1: "#D9EAF4",
    2: "#B8D8E8",
    3: "#8FC1DA",
    4: "#5FA6C8",
    5: "#347FA8",
}

# Collapse related subclasses in the same way as the Figure 2 class plot.
class_groups = {
    "Cer_AP": "Cer",
    "Cer_HDS": "Cer",
    "Cer_NDS": "Cer",
    "Cer_NP": "Cer",
    "Cer_NS": "Cer",
    "Cer_EODS": "Cer",
    "Cer_EOS": "Cer",
    "Cer_HS": "Cer",
    "Cer_AS": "Cer",
    "Cer_ADS": "Cer",
    "HexCer_AP": "HexCer",
    "HexCer_NDS": "HexCer",
    "HexCer_HS": "HexCer",
    "HexCer_NS": "HexCer",
    "HexCer_EOS": "HexCer",
    "HexCer_HDS": "HexCer",
    "AHexCer": "HexCer",
    "Hex2Cer": "HexCer",
    "Hex3Cer": "HexCer",
    "SHexCer": "HexCer",
    "ASM": "SM",
    "LNAPE": "PE",
    "LNAPS": "PS",
    "NAPE": "PE",
    "NAPS": "PS",
    "PE_Cer": "Cer",
    "PI_Cer": "Cer",

    "GM3": "GL",
    "DGDG": "DG",

    "MGDG": "DG",

}


def simplify_lipid_class(value):
    lipid_class = str(value).strip()
    lipid_class = class_groups.get(lipid_class, lipid_class)
    lipid_class = lipid_class.replace("Ox", "").replace("Ether", "")
    if lipid_class.startswith("L"):
        lipid_class = lipid_class[1:]
    return class_groups.get(lipid_class, lipid_class)


def prepare_coverage(frame):
    df = frame.copy()

    required_columns = {
        "Dataset",
        "Metabolite name",
        "Ontology",
        "Overall Match Status",
    }
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            "The input workbook is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    df["Dataset"] = (
        df["Dataset"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )
    df["Ontology"] = (
        df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    true_values = [True, "TRUE", "True", "true", 1, "1"]
    confirmed_df = df[
        df["Dataset"].isin(datasets)
        & df["Overall Match Status"].isin(true_values)
    ].copy()

    # Use sum-composition names so coverage is comparable with Figures 1 and 2.
    confirmed_df["Lipid Name"] = (
        confirmed_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )
    confirmed_df["Lipid Class"] = confirmed_df["Ontology"].apply(
        simplify_lipid_class
    )

    # Count a lipid once per dataset, then assign every unique class/name pair to
    # the number of datasets in which it was confirmed.
    confirmed_df = confirmed_df.drop_duplicates(
        subset=["Dataset", "Lipid Class", "Lipid Name"]
    )
    unique_lipids = (
        confirmed_df
        .groupby(["Lipid Class", "Lipid Name"])["Dataset"]
        .nunique()
        .rename("Coverage")
        .reset_index()
    )

    coverage_groups = list(range(1, len(datasets) + 1))
    coverage_counts = pd.crosstab(
        unique_lipids["Lipid Class"],
        unique_lipids["Coverage"],
    ).reindex(columns=coverage_groups, fill_value=0).sort_index()

    class_totals = coverage_counts.sum(axis=1)
    coverage_percentages = coverage_counts.div(class_totals, axis=0) * 100

    return coverage_counts, coverage_percentages


def plot_coverage(coverage_percentages, ax):
    coverage_groups = list(range(1, len(datasets) + 1))
    y_positions = np.arange(len(coverage_percentages.index))
    left = np.zeros(len(coverage_percentages.index))

    for coverage in coverage_groups:
        values = coverage_percentages[coverage].to_numpy()
        ax.barh(
            y_positions,
            values,
            left=left,
            height=0.72,
            color=coverage_colors[coverage],
            edgecolor="white",
            linewidth=0.7,
            label=f"{coverage} Dataset{'s' if coverage != 1 else ''}",
        )
        left += values

    ax.set_yticks(y_positions, labels=coverage_percentages.index)
    ax.margins(y=0.01)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_xlabel("Confirmed Lipids in Class (%)", labelpad=4)
    ax.set_ylabel("Lipid Class")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        title="Dataset Coverage",
        title_fontsize=9,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
        ncol=1,
    )

    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)
    return ax


DATASET_LABELS = {
    "OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT",
    "Stellar": "Stellar", "Stellarboth": "Stellar\u00b1",
}


def prepare_overlap(source):
    """Apply Figure 2B's filtering, deduplication, and class simplification."""
    from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import count_number_of_datasets
    from PublicScripts.Lipids.LipidFunctions import hg_simplifier, simplify_class_df

    frame = source.copy()
    for column in ["Dataset", "Ontology"]:
        frame[column] = (
            frame[column].astype(str).str.replace("\u00a0", " ", regex=False).str.strip()
        )
    frame = frame.loc[
        frame["Dataset"].isin(datasets)
        & frame["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    ].copy()
    frame["Metabolite name"] = (
        frame["Metabolite name"].astype(str).str.split("|", n=1).str[0].str.strip()
    )
    frame = frame.drop_duplicates(["Dataset", "Ontology", "Metabolite name"])
    frame["Dataset"] = frame["Dataset"].replace(DATASET_LABELS)
    frame = frame.drop_duplicates(["Dataset", "Metabolite name"]).copy()
    frame["Ontology"] = frame["Ontology"].replace(hg_simplifier)
    frame = simplify_class_df(frame)
    return count_number_of_datasets(frame)


def plot_class_comparison(source):
    _, percentages = prepare_coverage(source)
    if percentages.empty:
        raise ValueError("No confirmed lipids in the five manuscript datasets.")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    plot_coverage(percentages, ax)
    fig.tight_layout(pad=0.5)
    return fig



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(topdir) / input_file)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    source = (pd.read_csv(args.input) if args.input.suffix.casefold() == ".csv"
              else pd.read_excel(args.input, sheet_name="Sheet1"))
    counts, percentages = prepare_coverage(source)
    fig = plot_class_comparison(source)
    output_dir = args.output_dir or args.input.parent / "Final Figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        fig.savefig(
            output_dir / f"{output_basename}.{extension}",
            dpi=600, bbox_inches="tight", pad_inches=0.04,
        )
    print("\nConfirmed lipids in each class:")
    print(counts.sum(axis=1).to_string())
    print("\nCoverage percentages within each lipid class:")
    print(percentages.round(1).to_string())
    if args.no_show:
        plt.close(fig)
    else:
        plt.show(block=True)


if __name__ == "__main__":
    main()
