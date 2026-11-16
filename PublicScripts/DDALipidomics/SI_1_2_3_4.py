"""Plot RT versus m/z before and after lipid validation, for both polarities.

SI 1 and 2 show raw positive and negative MS-DIAL CSV rows, respectively, from
the five manuscript datasets. SI 3 and 4 show positive and negative IDs from
Combined_Fragment_Assignment_Colored.xlsx accepted by Overall Match Status,
excluding provisional entries. No validation or annotation filtering is applied
to the raw CSVs; only rows with unusable RT or m/z are excluded.
Each plotted point is one input row: names are not deduplicated. Class families
and colors follow LipidFunctions; unrecognized classes are displayed as Other.

Running this script saves four PNG (600 dpi)/PDF figure pairs beside the
workbook. Raw CSVs default to the workbook directory. Use --input PATH,
--raw-dir PATH, --output-dir PATH, or --no-show as needed.
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from PublicScripts.Lipids.LipidFunctions import class_color_map, hg_simplifier


INPUT_FILE = Path(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline\Combined_Fragment_Assignment_Colored.xlsx")
DATASETS = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
DATASET_LABELS = {"OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT",
                  "Stellar": "Stellar", "Stellarboth": "Stellar\u00b1"}
TRUE_VALUES = [True, "TRUE", "True", "true", 1, "1"]
RT_COLUMN = "Average Rt(min)"
MZ_COLUMN = "Average Mz"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def clean_text(series):
    return series.fillna("").astype(str).str.replace(
        "\u00a0", " ", regex=False
    ).str.strip()


def prepare_lipids(frame, raw=False):
    required = {"Dataset", "Polarity", "Ontology", RT_COLUMN, MZ_COLUMN}
    if not raw:
        required.add("Overall Match Status")
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))
    frame = frame.copy()
    frame["Dataset"] = clean_text(frame["Dataset"]).str.casefold().map(
        {name.casefold(): name for name in DATASETS}
    )
    frame["Polarity"] = clean_text(frame["Polarity"]).str.casefold()
    frame = frame.loc[
        frame["Dataset"].notna() & frame["Polarity"].isin(["positive", "negative"])
    ].copy()
    for column in [RT_COLUMN, MZ_COLUMN]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    usable = (
        np.isfinite(frame[[RT_COLUMN, MZ_COLUMN]]).all(axis=1)
        & frame[RT_COLUMN].ge(0) & frame[MZ_COLUMN].gt(0)
    )
    print(f"Excluded {int((~usable).sum())} rows with unusable RT or m/z.")
    frame = frame.loc[usable].copy()
    if frame.empty:
        raise ValueError("No usable lipid IDs remain for the five selected datasets.")
    frame["Accepted"] = False if raw else frame["Overall Match Status"].isin(TRUE_VALUES)
    if not raw and "Validation Status" in frame:
        frame["Accepted"] &= clean_text(frame["Validation Status"]).str.casefold().ne("provisional")
    families = {name.casefold(): group for name, group in hg_simplifier.items()}
    known_classes = {name.casefold(): name for name in class_color_map}
    frame["Plot Class"] = clean_text(frame["Ontology"]).map(
        lambda name: families.get(name.casefold(), name)
    ).map(lambda name: known_classes.get(name.casefold(), "Other"))
    return frame


def load_raw_lipids(directory):
    frames = []
    for dataset in DATASETS:
        for suffix, polarity in [("pos", "positive"), ("neg", "negative")]:
            frame = pd.read_csv(directory / f"{dataset}_{suffix}.csv")
            frame["Dataset"] = dataset
            frame["Polarity"] = polarity
            frames.append(frame)
    return prepare_lipids(pd.concat(frames, ignore_index=True), raw=True)


def plot_scatter_panels(frame, polarity, filtered, rt_limit, mz_limit):
    mode_frame = frame.loc[frame["Polarity"].eq(polarity)]
    plotted = mode_frame.loc[mode_frame["Accepted"]] if filtered else mode_frame
    # The legend must describe the plotted population, including acceptance filtering.
    classes = sorted(set(plotted["Plot Class"]).difference({"Other"}))
    if plotted["Plot Class"].eq("Other").any():
        classes.append("Other")
    fig, axes = plt.subplots(
        3, 2, figsize=(6.5, 6.5), sharex=True, sharey=True, layout="constrained",
    )
    for label, dataset, ax in zip("ABCDE", DATASETS, axes.flat):
        group = plotted.loc[plotted["Dataset"].eq(dataset)]
        colors = group["Plot Class"].map(class_color_map).fillna("#777777")
        ax.scatter(group[RT_COLUMN], group[MZ_COLUMN], c=colors,
                   s=9, alpha=0.78, linewidths=0, rasterized=True)
        ax.set_title(DATASET_LABELS[dataset], pad=5)
        ax.set_xlim(0, rt_limit)
        ax.set_ylim(0, mz_limit)
        ax.set_xticks(np.arange(0, rt_limit + 0.1, 5))
        ax.set_yticks(np.arange(0, mz_limit + 0.1, 400))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1)
        ax.spines["bottom"].set_linewidth(1)
        ax.text(0.02, 1.04, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", ha="left", va="bottom")
        ax.text(0.96, 0.94, f"n = {len(group)}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)
    # The unused sixth cell holds the class key, matching the original figure.
    legend_ax = axes[2, 1]
    legend_ax.axis("off")
    if classes:
        legend_ax.legend(
            handles=[Line2D([], [], marker="o", linestyle="none",
                            color=class_color_map[name], markersize=4, label=name)
                     for name in classes],
            title="Lipid Class", title_fontsize=9, loc="center",
            ncol=3, frameon=False, handlelength=0.7, handletextpad=0.4,
            columnspacing=0.7, labelspacing=0.15,
        )
    for ax in axes[:, 0]:
        ax.set_ylabel("m/z", fontsize=9, fontstyle="italic", labelpad=3)
    # Label the lowest data axis in each column; the right column ends at D.
    for ax in [axes[2, 0], axes[1, 1]]:
        ax.set_xlabel("Retention Time (min)", fontsize=9, labelpad=3)
        ax.tick_params(axis="x", labelbottom=True)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--raw-dir", type=Path,
                        help="Directory with the original DATASET_pos/neg.csv exports")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    frame = prepare_lipids(pd.read_excel(args.input, sheet_name="Sheet1"))
    raw_frame = load_raw_lipids(args.raw_dir or args.input.parent)
    output_dir = args.output_dir or args.input.parent / "Final Figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    # Use shared limits for all four plots and include high-mass IDs that were
    # clipped by the original figure's fixed 1600 m/z limit.
    all_rows = pd.concat([raw_frame, frame], ignore_index=True)
    rt_limit = max(20, int(np.ceil(all_rows[RT_COLUMN].max() / 5)) * 5)
    mz_limit = max(1600, int(np.ceil(all_rows[MZ_COLUMN].max() / 200)) * 200)
    counts = []
    for filtered in [False, True]:
        source = frame if filtered else raw_frame
        for polarity in ["positive", "negative"]:
            figure_number = 1 + 2 * int(filtered) + int(polarity == "negative")
            fig = plot_scatter_panels(source, polarity, filtered, rt_limit, mz_limit)
            stage = "Post_Filtered" if filtered else "Pre_Filtered"
            basename = f"SI_{figure_number}"
            for extension in ["png", "pdf"]:
                path = output_dir / f"{basename}.{extension}"
                fig.savefig(path, dpi=600)
                print(f"Saved {path}")
            selected = source.loc[source["Polarity"].eq(polarity)]
            if filtered:
                selected = selected.loc[selected["Accepted"]]
            for dataset in DATASETS:
                counts.append({"Stage": stage, "Polarity": polarity,
                               "Dataset": dataset,
                               "IDs": int(selected["Dataset"].eq(dataset).sum())})
            if args.no_show:
                plt.close(fig)
    summary = pd.DataFrame(counts)
    print("\nPlotted ID counts:")
    print(summary.to_string(index=False))
    if not args.no_show:
        plt.show(block=True)


if __name__ == "__main__":
    main()
