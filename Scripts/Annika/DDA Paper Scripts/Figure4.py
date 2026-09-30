"""Plot dataset coverage of lipids containing each tail, using Figure 3's format.

Default: confirmed IDs from the five manuscript datasets. Each distinct class
and chain composition counts once per dataset across features/adducts/polarities.
Each tail has a 100% stacked horizontal bar: segments show the fraction of
unique lipids containing that tail found in exactly 1, 2, 3, 4, or 5 datasets.
A lipid containing two copies of a tail counts once in that tail's bar.
Coverage uses resolved chain compositions so tail identities remain meaningful.

Show tails reaching 1% of all annotated tail copies
in at least one selected dataset, ordered by carbon number and double-bond
count. This frequency cutoff only selects the tails to show; bar percentages
describe dataset coverage. Other tails, LCB entries, and O- ether tails are hidden.

An additional figure uses Figure 2's stacked-bar format: lipids shared across
5, 4, 3, or 2 datasets, followed by lipids unique to each instrument/mode.
Colored segments count annotations for each displayed tail. Each unique lipid
contributes once per distinct tail, so stack totals count lipid-tail pairs.
Rarity and instrument detection tables are also exported for all eligible tails.
Rarity means occurrence among unique annotated lipids, not measured MS1 signal.

Sum-only multi-chain annotations are excluded, rather than interpreting a total
such as PC 34:1 as one 34:1 tail. Chain modifications and sphingoid backbones
remain distinct. These are annotation compositions, not abundance measurements
or independent confirmation of every chain by fragment matching.

Run with --min-tail-percent 1 (default) or --no-show. --tail 18:0 selects the
tail coverage printed in the console; it does not override the figure's threshold. Use --all-ids
explicitly to analyze an unvalidated CSV such as Combined_DDA.csv from DDA prep.
"""

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np
import pandas as pd


TOPDIR = Path(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline")
DATASETS = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
DATASET_LABELS = {
    "OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT",
    "Stellar": "Stellar", "Stellarboth": "Stellar (+/-)",
}
RARITY_BINS = [0, 0.5, 1, 2, 5, np.inf]
RARITY_LABELS = ["\u22640.5%", "0.5\u20131%", "1\u20132%", "2\u20135%", ">5%"]
# Same sequential blue palette as Figure 3.
COVERAGE_COLORS = {
    1: "#D9EAF4",
    2: "#B8D8E8",
    3: "#8FC1DA",
    4: "#5FA6C8",
    5: "#347FA8",
}

# Expected chains in the core annotation. Extra (FA ...) or leading (O-...)
# substituents are parsed separately. Unknown classes are reported, not guessed.
CHAIN_COUNTS = {
    **dict.fromkeys(["FA", "CE", "CAR", "NAE", "MG", "LPC", "LPE", "LPG", "LPI", "LPS", "LPA"], 1),
    **dict.fromkeys(["PC", "PE", "PG", "PI", "PS", "PA", "DG", "Cer", "SM", "HexCer", "Hex2Cer", "Hex3Cer", "GM3", "DGDG", "MGDG", "DGGA", "SQDG", "DLCL"], 2),
    **dict.fromkeys(["TG", "ADGGA", "HBMP", "MLCL"], 3),
    "CL": 4,
}
SPHINGOLIPIDS = {"Cer", "SM", "HexCer", "Hex2Cer", "Hex3Cer", "GM3"}
CHAIN_PATTERN = re.compile(r"(?P<prefix>[OP]-)?(?P<carbon>\d+):(?P<db>\d+)(?P<mod>;O\d*|\(\d+OH\))?")

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def parse_lipid_tails(name, lipid_class):
    """Return (detailed class, chain tuple, exclusion reason).

    Parse the current name, rather than T1-T4: those columns can be stale after
    a name edit, and T1 also contains unresolved totals in the pipeline output.
    """
    if not isinstance(name, str) or not name.strip():
        return "", (), "Missing lipid name"
    base_class = str(lipid_class).strip().removeprefix("Ether").removeprefix("Ox")
    if base_class not in CHAIN_COUNTS:
        return "", (), "Unsupported class or no discrete lipid tails"
    detailed = name.split("|")[-1].strip().split(maxsplit=1)
    if len(detailed) != 2:
        return "", (), "Missing chain annotation"
    detailed_class, body = detailed
    extra = []

    def collect_extra(match):
        extra.append(match.group(1).strip())
        return ""

    body = re.sub(r"^\(([^()]+)\)", collect_extra, body)
    body = re.sub(r"\(FA\s+([^()]+)\)", collect_extra, body)
    core = [part.strip() for part in re.split(r"[_/]", body)]
    core = [part for part in core if part != "0:0"]  # Lyso empty-chain placeholder.
    if len(core) != CHAIN_COUNTS[base_class]:
        return detailed_class, (), "Unresolved total or incomplete chain annotation"
    chains = []
    for position, chain in enumerate(core + extra):
        match = CHAIN_PATTERN.fullmatch(chain)
        if match is None or int(match['carbon']) == 0:
            return detailed_class, (), "Unsupported chain notation"
        value = f"{match['prefix'] or ''}{int(match['carbon'])}:{int(match['db'])}{match['mod'] or ''}"
        if base_class in SPHINGOLIPIDS and position == 0:
            value = "LCB " + value
        chains.append(value)
    # Tail identity comparison ignores sn-position order; retain multiplicities
    # in the species key so TG 16:0_18:0_18:0 and 16:0_16:0_18:0 stay distinct.
    return detailed_class, tuple(sorted(chains)), ""


def prepare_lipids(input_path, datasets, all_ids=False):
    source = pd.read_csv(input_path) if input_path.suffix.lower() == ".csv" else pd.read_excel(input_path)
    required = {"Dataset", "Ontology", "Metabolite name"}
    if not all_ids:
        required.add("Overall Match Status")
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Missing columns: " + ", ".join(sorted(missing)) + ". Raw DDA annotations require --all-ids.")
    for column in ["Dataset", "Ontology", "Metabolite name"]:
        source[column] = source[column].fillna("").astype(str).str.replace("\u00a0", " ", regex=False).str.strip()
    source["Dataset"] = source["Dataset"].str.casefold().map({d.casefold(): d for d in datasets})
    source = source.loc[source["Dataset"].notna()].copy()
    if not all_ids:
        source = source.loc[source["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])].copy()
    if source.empty:
        raise ValueError("No eligible IDs in the selected datasets.")
    audit = source[["Dataset", "Ontology", "Metabolite name"]].copy()
    parsed = [parse_lipid_tails(row["Metabolite name"], row["Ontology"]) for _, row in source.iterrows()]
    audit["Detailed class"] = [item[0] for item in parsed]
    audit["Chains"] = [item[1] for item in parsed]
    audit["Exclusion reason"] = [item[2] for item in parsed]
    audit["Species"] = audit["Detailed class"] + " " + audit["Chains"].map("_".join)
    resolved = audit["Exclusion reason"].eq("")
    duplicate_indices = audit.loc[resolved].loc[
        audit.loc[resolved].duplicated(["Dataset", "Ontology", "Species"])
    ].index
    audit.loc[duplicate_indices, "Exclusion reason"] = "Duplicate chain composition within dataset"
    lipids = audit.loc[audit["Exclusion reason"].eq("")].copy()
    if lipids.empty:
        raise ValueError("No resolved chain annotations were found in the selected IDs.")
    lipids = lipids.reset_index(drop=True)
    occurrences = pd.DataFrame([
        {"Dataset": row["Dataset"], "Ontology": row["Ontology"], "Species": row["Species"],
         "Tail": tail, "Tail copies": row["Chains"].count(tail)}
        for _, row in lipids.iterrows() for tail in sorted(set(row["Chains"]))
    ])
    return lipids, occurrences, audit


def frequency_summary(lipids, occurrences, group_columns):
    """Return lipid prevalence and tail-copy composition, including absent tails."""
    denominators = lipids.groupby(group_columns).size().rename("Resolved lipids")
    counts = occurrences.groupby(group_columns + ["Tail"]).size().unstack("Tail", fill_value=0)
    counts = counts.reindex(denominators.index, fill_value=0)
    summary = counts.stack().rename("Lipids with tail").reset_index()
    summary = summary.merge(denominators.reset_index(), on=group_columns, validate="many_to_one")
    copies = occurrences.groupby(group_columns + ["Tail"])["Tail copies"].sum().rename("Tail occurrences")
    summary = summary.merge(copies.reset_index(), on=group_columns + ["Tail"], how="left", validate="one_to_one")
    summary["Tail occurrences"] = summary["Tail occurrences"].fillna(0).astype(int)
    summary["Lipids with tail (%)"] = 100 * summary["Lipids with tail"] / summary["Resolved lipids"]
    summary["Total tail occurrences"] = summary.groupby(group_columns)["Tail occurrences"].transform("sum")
    summary["Tail composition (%)"] = (
        100 * summary["Tail occurrences"] / summary["Total tail occurrences"].replace(0, np.nan)
    )
    return summary


def save_figure(fig, output_dir, name):
    for extension in ["png", "pdf"]:
        path = output_dir / f"{name}.{extension}"
        fig.savefig(path, dpi=600)
        print(f"Saved {path}")


def tail_structural_type(tail):
    """Order fatty-acyl chains, ether-linked chains, LCBs, then the remainder."""
    if tail == "Other tails":
        return 3
    if tail.startswith("LCB "):
        return 2
    if tail.startswith(("O-", "P-")):
        return 1
    return 0


def tail_sort_key(tail):
    """Sort by structural type, carbon count, and unsaturation."""
    structural_type = tail_structural_type(tail)
    match = CHAIN_PATTERN.fullmatch(tail.removeprefix("LCB "))
    if match is None:
        return (structural_type, np.inf, np.inf, tail)
    return (structural_type, int(match["carbon"]), int(match["db"]), tail)


def tail_composition_table(summary, datasets, min_tail_percent=1.0):
    """Hide Other, LCB, and O- tails without changing the full-distribution percentages."""
    counts = summary.pivot(index="Tail", columns="Dataset", values="Tail occurrences")
    counts = counts.reindex(columns=datasets).fillna(0)
    denominators = counts.sum(axis=0)
    frequencies = counts.div(denominators.replace(0, np.nan), axis=1) * 100
    eligible = (
        frequencies.max(axis=1).ge(min_tail_percent)
        & ~frequencies.index.str.startswith(("LCB ", "O-"))
        & (frequencies.index != "Other tails")
    )
    tails = sorted(frequencies.index[eligible], key=tail_sort_key)
    return frequencies.reindex(index=tails)


def tail_coverage_tables(occurrences, datasets, tails):
    """Count unique lipids per tail by the number of selected datasets containing them."""
    # Ignore tail-copy multiplicity, adducts, and duplicate features when counting
    # dataset coverage. Each class/chain composition contributes once per tail.
    selected = occurrences.loc[
        occurrences["Dataset"].isin(datasets) & occurrences["Tail"].isin(tails)
    ]
    lipid_coverage = (
        selected.groupby(["Ontology", "Species", "Tail"])["Dataset"]
        .nunique().rename("Dataset Coverage").reset_index()
    )
    coverage_groups = list(range(1, len(datasets) + 1))
    counts = pd.crosstab(lipid_coverage["Tail"], lipid_coverage["Dataset Coverage"])
    counts = counts.reindex(index=tails, columns=coverage_groups, fill_value=0)
    totals = counts.sum(axis=1)
    percentages = counts.div(totals.replace(0, np.nan), axis=0) * 100
    return counts, percentages, lipid_coverage


def plot_coverage(percentages, datasets, population, ax):
    """Use Figure 3's 100% horizontal stacks and dataset-coverage legend."""
    y_positions = np.arange(len(percentages))
    left = np.zeros(len(percentages))
    colors = list(COVERAGE_COLORS.values())
    if len(datasets) > len(colors):
        colors = plt.get_cmap("Blues")(np.linspace(0.15, 0.75, len(datasets)))
    for coverage in percentages.columns:
        values = percentages[coverage].to_numpy()
        ax.barh(
            y_positions, values, left=left, height=0.72,
            color=colors[coverage - 1], edgecolor="white", linewidth=0.7,
            label=f"{coverage} Dataset{'s' if coverage != 1 else ''}",
        )
        left += values
    ax.set_yticks(y_positions, labels=percentages.index)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 101, 20))
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.tick_params(axis="both", labelsize=8)
    ax.set_xlabel(f"{population} Lipids Containing Tail (%)", labelpad=12, fontsize=9)
    ax.set_ylabel("Tail Identity", fontsize=9)
    ax.set_title(f"Dataset Coverage of {population} Lipids by Tail", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    for row in range(1, len(percentages)):
        if tail_structural_type(percentages.index[row]) != tail_structural_type(percentages.index[row - 1]):
            ax.axhline(row - 0.5, color="#999999", linewidth=0.7)
    if percentages.empty:
        ax.text(0.5, 0.5, "No eligible tails meet the threshold",
                transform=ax.transAxes, ha="center", va="center", fontsize=8, wrap=True)
    ax.legend(
        title="Dataset Coverage", title_fontsize=9, fontsize=8,
        bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, ncol=1,
    )
    return ax


def instrument_tail_tables(occurrences, datasets):
    """Compare detection against one pooled lipid catalog shared by all datasets."""
    selected = occurrences.loc[occurrences["Dataset"].isin(datasets)].drop_duplicates(
        ["Dataset", "Ontology", "Species", "Tail"]
    )
    present_datasets = set(selected["Dataset"])
    total_lipids = len(selected[["Ontology", "Species"]].drop_duplicates())
    selected = selected.loc[
        ~selected["Tail"].str.startswith(("LCB ", "O-")) & selected["Tail"].ne("Other tails")
    ]
    pooled = selected.drop_duplicates(["Ontology", "Species", "Tail"]).groupby("Tail").size()
    tails = sorted(pooled.index, key=tail_sort_key)
    counts = pd.crosstab(selected["Tail"], selected["Dataset"]).reindex(
        index=tails, columns=datasets, fill_value=0
    )
    # A selected dataset with no resolved lipids has no detection data.
    for dataset in datasets:
        if dataset not in present_datasets:
            counts[dataset] = np.nan
    metadata = pooled.reindex(tails).rename("Pooled lipids with tail").to_frame()
    metadata["Total pooled lipids"] = total_lipids
    metadata["Pooled lipid prevalence (%)"] = (
        100 * metadata["Pooled lipids with tail"] / (total_lipids or np.nan)
    )
    metadata["Prevalence bin"] = pd.cut(
        metadata["Pooled lipid prevalence (%)"], RARITY_BINS, labels=RARITY_LABELS
    )
    percentages = counts.div(metadata["Pooled lipids with tail"], axis=0) * 100
    return counts, percentages, metadata


def rarity_detection_tables(counts, metadata):
    """Count distinct tails detected at least once, within pooled rarity bins."""
    bin_totals = metadata.groupby("Prevalence bin", observed=False).size().rename("Pooled tail identities")
    detected = counts.gt(0).groupby(metadata["Prevalence bin"], observed=False).sum()
    for dataset in counts.columns:
        if counts[dataset].isna().all():
            detected[dataset] = np.nan
    percentages = detected.div(bin_totals.replace(0, np.nan), axis=0) * 100
    return detected, percentages, bin_totals


def tail_overlap_tables(occurrences, datasets, tails):
    """Assign each lipid-tail pair to one shared or instrument-only category."""
    selected = occurrences.loc[
        occurrences["Dataset"].isin(datasets) & occurrences["Tail"].isin(tails)
    ].drop_duplicates(["Dataset", "Ontology", "Species", "Tail"])
    membership = selected.groupby(["Ontology", "Species", "Tail"])["Dataset"].agg(
        **{"Dataset Coverage": "nunique", "Only Dataset": "first"}
    ).reset_index()
    shared = membership["Dataset Coverage"].gt(1)
    membership["Overlap Category"] = membership["Only Dataset"].map(
        lambda dataset: DATASET_LABELS.get(dataset, dataset) + " Only"
    )
    membership.loc[shared, "Overlap Category"] = (
        membership.loc[shared, "Dataset Coverage"].astype(str) + " Datasets"
    )
    membership.loc[shared, "Only Dataset"] = ""
    categories = [f"{coverage} Datasets" for coverage in range(len(datasets), 1, -1)]
    categories += [DATASET_LABELS.get(dataset, dataset) + " Only" for dataset in datasets]
    counts = pd.crosstab(membership["Overlap Category"], membership["Tail"]).reindex(
        index=categories, columns=tails, fill_value=0
    )
    return counts, membership


def plot_instrument_comparison(counts, population, ax):
    """Match Figure 2: stacked counts by overlap category, colored by tail."""
    x_positions = np.arange(len(counts))
    bottom = np.zeros(len(counts))
    palette = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    for index, tail in enumerate(counts.columns):
        values = counts[tail].to_numpy()
        ax.bar(x_positions, values, bottom=bottom, width=0.65,
               color=palette[index % len(palette)], edgecolor="white", linewidth=0.3,
               label=tail)
        bottom += values
    ax.set_xticks(x_positions, labels=counts.index, rotation=45, ha="right")
    ax.tick_params(axis="both", labelsize=8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(1, bottom.max() * 1.08) if len(bottom) else 1)
    ax.set_ylabel("Tail Annotations (count)", fontsize=9)
    ax.set_title(f"Tail Distribution by Overlap Category\n{population} Lipids", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    if len(counts.columns):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title="Tail Identity", title_fontsize=9,
                  fontsize=8, bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False,
                  ncol=2 if len(counts.columns) > 12 else 1)
    else:
        ax.text(0.5, 0.5, "No eligible tails meet the threshold", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, wrap=True)
    return ax


def plot_tail_comparison(percentages, counts, datasets, output_dir, population):
    """Stack the dataset-coverage and overlap-category panels vertically."""
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.5),
                             gridspec_kw={"height_ratios": [1.5, 1]})
    plot_coverage(percentages, datasets, population, axes[0])
    plot_instrument_comparison(counts, population, axes[1])
    for label, ax in zip(["A", "B"], axes):
        ax.text(-0.12, 1.04, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")
    fig.tight_layout(h_pad=2.5)
    save_figure(fig, output_dir, "Figure4")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=TOPDIR / "Combined_Fragment_Assignment_Colored.xlsx")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--tail", default="18:0")
    parser.add_argument(
        "--min-tail-percent", type=float, default=1.0,
        help="Show tails excluding LCB and O- entries that reach this percentage in any selected dataset (default: 1).",
    )
    parser.add_argument("--all-ids", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.min_tail_percent <= 100:
        parser.error("--min-tail-percent must be between 0 and 100")
    datasets = list(dict.fromkeys(args.datasets))
    output_dir = args.output_dir or args.input.parent
    lipids, occurrences, audit = prepare_lipids(args.input, datasets, args.all_ids)
    output_dir.mkdir(parents=True, exist_ok=True)
    audit.assign(Chains=audit.Chains.map("; ".join)).to_csv(output_dir / "Tail_Annotation_Audit.csv", index=False)
    lipids.assign(Chains=lipids.Chains.map("; ".join)).to_csv(output_dir / "Tail_Resolved_Lipids.csv", index=False)
    occurrences.to_csv(output_dir / "Tail_Occurrences.csv", index=False)
    overall = frequency_summary(lipids, occurrences, ["Dataset"])
    by_class = frequency_summary(lipids, occurrences, ["Ontology", "Dataset"])
    overall.to_csv(output_dir / "Tail_Frequencies_by_Dataset.csv", index=False)
    by_class.to_csv(output_dir / "Tail_Frequencies_by_Class.csv", index=False)
    population = "Annotated" if args.all_ids else "Accepted"
    selected_frequencies = tail_composition_table(overall, datasets, args.min_tail_percent)
    selected_frequencies.to_csv(output_dir / "Tail_Selection_Frequencies.csv", index_label="Tail")
    coverage_counts, coverage_percentages, lipid_coverage = tail_coverage_tables(
        occurrences, datasets, selected_frequencies.index
    )
    coverage_counts.to_csv(output_dir / "Tail_Coverage_Counts.csv", index_label="Tail")
    coverage_percentages.to_csv(output_dir / "Tail_Coverage_Percentages.csv", index_label="Tail")
    lipid_coverage.to_csv(output_dir / "Tail_Lipid_Dataset_Coverage.csv", index=False)
    instrument_counts, instrument_percentages, tail_metadata = instrument_tail_tables(occurrences, datasets)
    instrument_counts.to_csv(output_dir / "Tail_Instrument_Detection_Counts.csv", index_label="Tail")
    instrument_percentages.to_csv(output_dir / "Tail_Instrument_Detection_Percentages.csv", index_label="Tail")
    tail_metadata.to_csv(output_dir / "Tail_Instrument_Pooled_Lipids.csv", index_label="Tail")
    overlap_counts, overlap_membership = tail_overlap_tables(
        occurrences, datasets, selected_frequencies.index
    )
    overlap_counts.to_csv(output_dir / "Tail_Overlap_Category_Counts.csv")
    overlap_membership.to_csv(output_dir / "Tail_Overlap_Lipid_Membership.csv", index=False)
    plot_tail_comparison(
        coverage_percentages, overlap_counts, datasets, output_dir, population
    )
    rarity_counts, rarity_percentages, rarity_totals = rarity_detection_tables(instrument_counts, tail_metadata)
    rarity_counts.to_csv(output_dir / "Tail_Rarity_Detection_Counts.csv")
    rarity_percentages.to_csv(output_dir / "Tail_Rarity_Detection_Percentages.csv")
    rarity_totals.to_csv(output_dir / "Tail_Rarity_Pooled_Counts.csv")
    print("\nTail identities detected by instrument within each pooled-prevalence bin (%):")
    print(rarity_percentages.rename(
        index=lambda label: str(label).replace("\u2264", "<=").replace("\u2013", "-")
    ).round(1).to_string())
    print("\nAnnotation accounting:")
    print(audit["Exclusion reason"].replace("", "Included unique chain composition").value_counts().to_string())
    print(f"\nDataset coverage of lipids containing {args.tail}:")
    tail_counts, tail_percentages, _ = tail_coverage_tables(occurrences, datasets, [args.tail])
    print(pd.DataFrame({
        "Unique lipids": tail_counts.loc[args.tail],
        "Lipids (%)": tail_percentages.loc[args.tail],
    }).round(2).to_string())
    if args.no_show:
        plt.close("all")
    else:
        plt.show(block=True)


if __name__ == "__main__":
    main()
