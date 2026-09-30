"""Exploratory rejection-by-class and intensity-versus-coverage figures.

Run directly to save PNG/PDF figures and their CSV data beside the manuscript
workbook. Use --no-show to generate the files without opening plot windows.

The rejection figure uses the LipidChecker output and Figure 1's counting:
each class/sum-composition name counts once per dataset within each match-status
group. A name with both confirmed and rejected entries counts in both groups.
Rejection is defined by Overall Match Status, including failures without rules.
Separate removal workbooks and outlier predictions are not part of this count.
Intensity coverage uses confirmed sum-composition lipids. Each intensity point
is the median across lipids with that coverage, within one dataset/polarity.
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, PercentFormatter
import numpy as np
import pandas as pd


# Share Figure 2's lipid-class palette when this script is run directly.
repository_root = Path(__file__).resolve().parents[3]
if str(repository_root) not in sys.path:
    sys.path.insert(0, str(repository_root))

from PublicScripts.Lipids.LipidFunctions import class_color_map  # noqa: E402


TOPDIR = Path(
    r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
)
INPUT_FILE = "Combined_Fragment_Assignment_Colored.xlsx"
DATASETS = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
DATASET_LABELS = ["OT/OT", "OT/IT", "IT/IT", "Stellar", "Stellar (+/-)"]
DATASET_COLORS = ["#440154", "#414487", "#2A788E", "#22A884", "#7AD151"]
MIN_CLASS_IDS = 20  # Minimum confirmed + rejected IDs; all classes go to CSV.
INTENSITY_CLASSES = ["PC", "PE", "SM", "Cer", "HexCer", "TG"]

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def clean_text(series):
    return series.fillna("").astype(str).str.replace(
        "\u00a0", " ", regex=False
    ).str.strip()


def load_ids(input_path):
    frame = pd.read_excel(input_path, sheet_name="Sheet1")
    required = {
        "Dataset", "Ontology", "Metabolite name", "Polarity",
        "Overall Match Status", "Matched Required Fragments",
        "MS1 isotopic spectrum",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))

    for column in ["Dataset", "Ontology", "Metabolite name", "Polarity"]:
        frame[column] = clean_text(frame[column])
    frame["Dataset"] = frame["Dataset"].str.casefold().map(
        {name.casefold(): name for name in DATASETS}
    )
    frame = frame.loc[frame["Dataset"].notna()].copy()
    if frame.empty:
        raise ValueError("No ID entries found for the selected datasets.")
    frame["Ontology"] = frame["Ontology"].replace("", "Unknown")
    frame["Polarity"] = frame["Polarity"].str.capitalize()
    frame["Lipid Name"] = frame["Metabolite name"].str.split(
        "|", n=1
    ).str[0].str.strip()

    # Match Figure1.py exactly: its rejected population is the complement of
    # this accepted-value list, with no additional validation-rule filter.
    true_values = [True, "TRUE", "True", "true", 1, "1"]
    frame["Confirmed"] = frame["Overall Match Status"].isin(true_values)
    frame["Rejected"] = ~frame["Confirmed"]
    # LipidChecker.check_fragments returns False when no class/adduct rule exists.
    frame["No Rule"] = clean_text(
        frame["Matched Required Fragments"]
    ).str.casefold().str.startswith("no required fragments found for this type")
    return frame


def rejection_summary(frame, keys):
    """Count each lipid within dataset/status, exactly as Figure 1 does.

    Include status in the deduplication key so a name that has both passing
    and failing entries is represented once in each population. Using any()
    confirmation across a name would silently discard its rejected IDs.
    """
    ids = frame.drop_duplicates(["Dataset", "Ontology", "Lipid Name", "Confirmed"])
    summary = ids.groupby(keys).agg(**{
        "Total IDs": ("Rejected", "size"),
        "Confirmed IDs": ("Confirmed", "sum"),
        "Rejected IDs": ("Rejected", "sum"),
    })
    summary["Rejected (%)"] = 100 * summary["Rejected IDs"] / summary["Total IDs"]
    # Preserve rule availability as source-row metadata, not another filter.
    # A name may have multiple rows with different rule availability.
    source_counts = frame.groupby(keys).agg(**{
        "Source rows": ("Rejected", "size"),
        "Source rows without rules": ("No Rule", "sum"),
    })
    return summary.join(source_counts)


def save_figure(fig, output_dir, basename):
    for extension in ["png", "pdf"]:
        path = output_dir / f"{basename}.{extension}"
        fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.1)
        print(f"Saved {path}")


def plot_rejection_by_class(frame, output_dir, minimum_ids=MIN_CLASS_IDS):
    summary = rejection_summary(frame, ["Ontology"]).sort_values(
        ["Rejected (%)", "Total IDs"],
        ascending=[False, False], kind="stable",
    )
    by_dataset = rejection_summary(frame, ["Ontology", "Dataset"])
    dataset_totals = rejection_summary(frame, ["Dataset"]).reindex(DATASETS)
    summary["Shown in figure"] = summary["Total IDs"].ge(minimum_ids)
    summary.to_csv(output_dir / "Rejection_by_Lipid_Class.csv")
    by_dataset.to_csv(output_dir / "Rejection_by_Lipid_Class_and_Dataset.csv")
    dataset_totals.to_csv(output_dir / "Rejection_by_Dataset.csv")
    shown = summary.loc[summary["Shown in figure"]]
    if shown.empty:
        print(f"No classes have at least {minimum_ids} IDs; saved CSVs only.")
        return summary

    fig, ax = plt.subplots(figsize=(6.5, max(4.5, 0.23 * len(shown) + 1.4)))
    positions = np.arange(len(shown))
    rates = shown["Rejected (%)"]
    bars = ax.barh(positions, rates, color="#0072B2", height=0.7)
    for bar, rate in zip(bars, rates):
        # Keep 100% labels inside the plotting area.
        ax.text(
            rate - 1.5 if rate >= 90 else rate + 1.5,
            bar.get_y() + bar.get_height() / 2,
            f"{rate:.1f}%", va="center",
            ha="right" if rate >= 90 else "left",
            color="white" if rate >= 90 else "black", fontsize=8,
        )
    no_rules = shown["Source rows without rules"].eq(shown["Source rows"])
    ax.set_yticks(
        positions,
        labels=[
            f"{name}{'†' if no_rules.loc[name] else ''}  (n={int(n)})"
            for name, n in shown["Total IDs"].items()
        ],
    )
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_xlabel("IDs rejected (%)")
    ax.set_ylabel("Lipid class")
    ax.set_title("LipidChecker Rejection by Class", pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.15)
    ax.set_axisbelow(True)
    fig.text(
        0.5, 0.015,
        "Percent rejected = rejected / (confirmed + rejected), using Figure 1's counts.\n"
        "Each class/name counts once per dataset and status; a name can occur in both status groups.\n"
        f"Classes with n >= {minimum_ids} shown. \u2020 All source rows lack rules; LipidChecker marks them rejected.",
        ha="center", va="bottom", fontsize=7,
    )
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save_figure(fig, output_dir, "Rejected_ID_Percentage_by_Lipid_Class")
    print("\nCounts by dataset (same definition as Figure 1):")
    print(dataset_totals[["Confirmed IDs", "Rejected IDs", "Rejected (%)"]].round(1).to_string())
    print("\nPlotted class rejection rates:")
    print(shown[["Confirmed IDs", "Rejected IDs", "Rejected (%)"]].round(1).to_string())
    zero_classes = shown.index[shown["Rejected IDs"].eq(0)]
    print("\nNo rejected IDs in this LipidChecker output for: " + ", ".join(zero_classes))
    return summary


def first_ms1_intensity(spectrum):
    """Intensity of the first (M) isotope peak in the MS-DIAL export.

    This is a spectral peak intensity, not integrated chromatographic area or
    a sum of isotope/fragment peaks. Never substitute another isotope for M.
    """
    if not isinstance(spectrum, str) or not spectrum.strip():
        return np.nan
    try:
        mz_text, intensity_text = spectrum.split()[0].split(":")
        mz, intensity = float(mz_text), float(intensity_text)
    except (ValueError, IndexError):
        return np.nan
    return intensity if np.isfinite([mz, intensity]).all() and mz > 0 and intensity > 0 else np.nan


def prepare_intensity_data(frame):
    confirmed = frame.loc[frame["Confirmed"] & frame["Lipid Name"].ne("")].copy()
    lipid_key = ["Ontology", "Lipid Name"]
    # Recompute after dataset selection and confirmation filtering; the workbook
    # Dataset Count also includes rejected entries and datasets outside this set.
    # A confirmation in either polarity counts as presence in that dataset.
    coverage = confirmed.groupby(lipid_key)["Dataset"].nunique().rename("Dataset Count")
    confirmed["MS1 Intensity"] = confirmed["MS1 isotopic spectrum"].map(first_ms1_intensity)
    invalid = confirmed["MS1 Intensity"].isna()
    unknown_polarity = ~confirmed["Polarity"].isin(["Positive", "Negative"])
    print(f"Intensity: excluded {invalid.sum()} confirmed entries with invalid/nonpositive MS1 intensity.")
    print(f"Intensity: excluded {unknown_polarity.sum()} confirmed entries with unknown polarity.")
    # Coverage is calculated before excluding missing intensities: missing signal
    # measurements must not change how many datasets confirmed a lipid.
    intensity = confirmed.loc[~invalid & ~unknown_polarity].groupby(
        ["Dataset", "Polarity"] + lipid_key
    ).agg(**{
        "MS1 Intensity": ("MS1 Intensity", "median"),
        "Contributing ID entries": ("MS1 Intensity", "size"),
    }).reset_index().merge(coverage.reset_index(), on=lipid_key, validate="many_to_one")
    intensity["log10 MS1 Intensity"] = np.log10(intensity["MS1 Intensity"])
    return intensity


def plot_intensity_by_coverage(frame, output_dir):
    intensity = prepare_intensity_data(frame)
    intensity.to_csv(output_dir / "Intensity_vs_Dataset_Count_Lipids.csv", index=False)
    if intensity.empty:
        print("No confirmed lipids with usable MS1 intensities; skipped intensity figure.")
        return intensity
    summary = intensity.groupby(["Dataset", "Polarity", "Dataset Count"]).agg(**{
        "Lipids": ("MS1 Intensity", "size"),
        "Median MS1 Intensity": ("MS1 Intensity", "median"),
        "Q1 MS1 Intensity": ("MS1 Intensity", lambda values: values.quantile(0.25)),
        "Q3 MS1 Intensity": ("MS1 Intensity", lambda values: values.quantile(0.75)),
    })
    summary.to_csv(output_dir / "Intensity_vs_Dataset_Count_Summary.csv")
    fig, axes = plt.subplots(ncols=2, figsize=(7.5, 4.2), sharey=True)
    coverage_groups = list(range(1, len(DATASETS) + 1))
    markers = ["o", "s", "^", "D", "v"]
    for ax, polarity in zip(axes, ["Positive", "Negative"]):
        for dataset, label, color, marker in zip(DATASETS, DATASET_LABELS, DATASET_COLORS, markers):
            if (dataset, polarity) not in summary.index.droplevel("Dataset Count"):
                continue
            group = summary.loc[(dataset, polarity)].reindex(coverage_groups)
            medians = group["Median MS1 Intensity"]
            ax.plot(
                coverage_groups, medians, color=color, linewidth=1.5, label=label,
                linestyle="--" if dataset == "OTOT" else "-",
                zorder=4 if dataset == "OTOT" else 2,
            )
            # Hollow points identify medians based on fewer than 10 lipids.
            # Missing groups remain NaN so the line has a gap.
            for count, median, n in zip(coverage_groups, medians, group["Lipids"]):
                if pd.notna(median):
                    ax.plot(
                        count, median, marker=marker, linestyle="none", color=color,
                        markerfacecolor=color if n >= 10 else "white", markersize=5,
                        zorder=5 if dataset == "OTOT" else 3,
                    )
            if dataset == "OTOT" and polarity == "Negative" and pd.isna(medians.loc[1]):
                ax.text(
                    0.02, 0.98, "OT/OT: no lipids in the 1-dataset group",
                    transform=ax.transAxes, va="top", fontsize=7, color=color,
                )
        ax.set_title(polarity + " mode")
        ax.set_xticks(coverage_groups)
        ax.set_xlim(0.8, len(DATASETS) + 0.2)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:,.0f}"))
        ax.minorticks_off()
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Median MS1 intensity (a.u.; log scale)")
    fig.suptitle("Intensity versus Number of Datasets", fontsize=12)
    fig.supxlabel("Number of datasets confirming the lipid", y=0.23, fontsize=9)
    fig.legend(
        handles=[
            Line2D(
                [0], [0], color=color, marker=marker, markersize=5, label=label,
                linestyle="--" if dataset == "OTOT" else "-",
            )
            for dataset, label, color, marker in zip(DATASETS, DATASET_LABELS, DATASET_COLORS, markers)
        ], loc="lower center", bbox_to_anchor=(0.5, 0.13), ncol=5, frameon=False,
    )
    fig.text(
        0.5, 0.015,
        "Each point is the median intensity of lipids confirmed in that many datasets.\n"
        "Intensity = first MS1 isotope-peak signal; each lipid contributes once per dataset/polarity.\n"
        "Confirmation in either polarity counts toward dataset number. Hollow points: fewer than 10 lipids.",
        ha="center", va="bottom", fontsize=7,
    )
    fig.tight_layout(rect=(0, 0.27, 1, 0.92))
    save_figure(fig, output_dir, "MS1_Intensity_vs_Dataset_Count")
    return intensity


def plot_intensity_by_class(intensity, output_dir, classes=INTENSITY_CLASSES):
    """One PDF page per dataset; colored curves compare classes within a mode."""
    if intensity.empty:
        return
    summary = intensity.groupby(
        ["Dataset", "Polarity", "Ontology", "Dataset Count"]
    ).agg(**{
        "Lipids": ("MS1 Intensity", "size"),
        "Median MS1 Intensity": ("MS1 Intensity", "median"),
    })
    summary.to_csv(output_dir / "Intensity_vs_Dataset_Count_by_Class.csv")
    coverage_groups = list(range(1, len(DATASETS) + 1))
    classes = list(dict.fromkeys(classes))
    colors = [class_color_map.get(lipid_class, "#333333") for lipid_class in classes]
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h"]
    pdf_path = output_dir / "MS1_Intensity_vs_Dataset_Count_by_Class.pdf"
    with PdfPages(pdf_path) as pdf:
        for dataset, dataset_label in zip(DATASETS, DATASET_LABELS):
            fig, axes = plt.subplots(ncols=2, figsize=(7.5, 4.2), sharey=True)
            for ax, polarity in zip(axes, ["Positive", "Negative"]):
                present = False
                for index, (lipid_class, color) in enumerate(zip(classes, colors)):
                    if (dataset, polarity, lipid_class) not in summary.index.droplevel("Dataset Count"):
                        continue
                    group = summary.loc[(dataset, polarity, lipid_class)].reindex(coverage_groups)
                    medians = group["Median MS1 Intensity"]
                    present = True
                    ax.plot(coverage_groups, medians, color=color, linewidth=1.3)
                    for count, median, n in zip(coverage_groups, medians, group["Lipids"]):
                        if pd.notna(median):
                            ax.plot(
                                count, median, marker=markers[index % len(markers)],
                                linestyle="none", color=color, markersize=5,
                                markerfacecolor=color if n >= 5 else "white",
                            )
                ax.set_title(polarity + " mode")
                ax.set_xticks(coverage_groups)
                ax.set_xlim(0.8, len(DATASETS) + 0.2)
                ax.set_yscale("log")
                ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:,.0f}"))
                ax.minorticks_off()
                ax.spines[["top", "right"]].set_visible(False)
                ax.grid(axis="y", alpha=0.15)
                if not present:
                    ax.text(0.5, 0.5, "No confirmed lipids\nfor these classes", transform=ax.transAxes, ha="center")
            axes[0].set_ylabel("Median MS1 peak intensity\nin this dataset (a.u.; log scale)")
            fig.suptitle(f"{dataset_label}: Intensity by Lipid Class", fontsize=12)
            fig.supxlabel(
                f"Datasets confirming each lipid (out of {len(DATASETS)})", y=0.23, fontsize=9
            )
            fig.legend(
                handles=[
                    Line2D([0], [0], color=color, marker=markers[index % len(markers)], label=lipid_class)
                    for index, (lipid_class, color) in enumerate(zip(classes, colors))
                ], loc="lower center", bbox_to_anchor=(0.5, 0.13), ncol=min(6, len(classes)), frameon=False,
            )
            fig.text(
                0.5, 0.015,
                "Each point: class median for lipids confirmed in exactly that many datasets.\n"
                "Intensity comes from the dataset and mode shown; confirmations in either mode count toward x.\n"
                "Different x values group different lipids. Hollow points: fewer than 5 lipids. Missing groups have no point.",
                ha="center", va="bottom", fontsize=7,
            )
            fig.tight_layout(rect=(0, 0.27, 1, 0.92))
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.1)
            fig.savefig(output_dir / f"MS1_Intensity_by_Class_{dataset}.png", dpi=600, bbox_inches="tight")
    print(f"Saved class comparisons: {pdf_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=TOPDIR / INPUT_FILE)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--min-class-ids", type=int, default=MIN_CLASS_IDS)
    parser.add_argument("--intensity-classes", nargs="+", default=INTENSITY_CLASSES)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if args.min_class_ids < 1:
        parser.error("--min-class-ids must be at least 1")
    output_dir = args.output_dir or args.input.parent
    frame = load_ids(args.input)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Included {len(frame)} ID entries from {frame['Dataset'].nunique()} selected datasets.")
    print(f"Source workbook: {args.input.resolve()}")
    print(f"Rejected source rows: {frame['Rejected'].sum()}; rows without rules: {frame['No Rule'].sum()}.")
    plot_rejection_by_class(frame, output_dir, args.min_class_ids)
    intensity = plot_intensity_by_coverage(frame, output_dir)
    plot_intensity_by_class(intensity, output_dir, args.intensity_classes)
    if args.no_show:
        plt.close("all")
    else:
        plt.show(block=True)


if __name__ == "__main__":
    main()
