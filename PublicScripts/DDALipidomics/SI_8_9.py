"""Compare adduct-specific class networks and within-adduct normalized strengths.

Rebuild networks from MS/MS spectrum and Average Mz, separately for each dataset,
ORIGINAL ontology, polarity and adduct BEFORE calculating strength. Never reuse
the input's pooled Class Strength. Use the existing network algorithm: mean of
fragment and neutral-loss cosine similarities, edges >0.1, weighted degree as
strength. Preserve each network's Tolerance (default 0.05 Da if absent). Networks
use accepted, non-provisional usable spectra by default; groups with fewer
than three distinct lipid identities are not scored. No cross-adduct edges exist.

Then select the plotted non-provisional identifications and
collapse duplicate features to their median within dataset/original-class/
polarity/adduct/lipid identity (name before '|', ignoring case and whitespace).
Normalize each median to the maximum within that SAME dataset/original-class/
polarity/adduct and plotted status population. Every nonzero normalization group
must reach 1; all-zero groups remain zero. --all-ids combines statuses before
network construction, consolidation and normalization. --excluded-only retains
the exploratory all-candidate network, then plots excluded IDs only.
Missing/invalid spectra are not scored.

Keep original ontology classes represented by accepted IDs, with separate
adduct panels. Never combine any original classes: PC/LPC/EtherPC,
PE/LPE/EtherPE, PG/LPG/EtherPG, PI/LPI/EtherPI, PS/LPS/EtherPS,
and PA/LPA/EtherPA all remain separate, as do other subclasses.
By default require the maximum eligible unique-lipid count across the five
datasets to be >=10, separately for each class/adduct/polarity and plotted status
population. Every dataset must have at least one eligible lipid; missing or
all-zero dataset groups exclude the entire panel. Exactly 10 qualifies. Override the
threshold with --min-lipids (>=2).
All eligible lipid identities contribute; no shared-lipid requirement.

For each retained class/adduct/polarity, estimate instrument means by OLS with
CR1 covariance clustered by lipid identity. Test all ten planned pairwise mean
contrasts using two-sided Student t, G-1 degrees of freedom (G lipid clusters).
Holm-adjust all planned comparisons across retained panels and both polarities,
separately per population, including unavailable comparisons as p=1. Export
pointwise 95% CIs. Significant pairs have brackets and adjusted-p stars.
Tests concern means; boxes show medians/IQRs with 1.5-IQR whiskers and all points.
P-values remain exploratory: clustering does not remove dependencies between
different lipids induced by the network or establish instrument-run replication.

Default: export accepted IDs only, up to twenty panels in four columns
per figure, PNG (600 dpi) and PDF. Figures are 6.5 inches wide and at most
7.5 inches tall, leaving 1.5 inches for a caption on a letter page with
one-inch margins. Shorter figures shrink to the occupied rows.
Panel letters continue across
pages within each polarity; positive and negative each start at A.
Recompute feature scores in memory
from Combined_Fragment_Assignment_Colored.xlsx; no intermediate score file is
required or saved. Save only PNG/PDF figures in Final Figures beside the input. --input accepts a CSV/workbook containing spectra, not just
precomputed scores. --output-dir overrides the destination; --no-show is headless.
"""

import argparse
from itertools import combinations
import os
from pathlib import Path
import sys

import matplotlib
if __name__ == "__main__":
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg" if "--no-show" in sys.argv else "WxAgg"))
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import t as t_distribution
from statsmodels.regression.linear_model import OLS

from Figure3 import datasets, topdir
from SI_Class_Heatmap import DATASET_LABELS

INPUT_FILE = Path(topdir) / "Combined_Fragment_Assignment_Colored.xlsx"
OUTPUT_BASENAMES = {"positive": "SI_8", "negative": "SI_9"}
NORMALIZED = "Within-Dataset Max-Normalized Class Strength"
POLARITY_COLORS = {"positive": "#347FA8", "negative": "#C43C39"}
FIGURE_SIZE = (6.5, 7.5)
MAX_PANELS_PER_PAGE = 20
MIN_LIPIDS = 10
NETWORK_GROUPS = ["Dataset", "Ontology", "Polarity", "Adduct type"]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 11, "xtick.labelsize": 8, "ytick.labelsize": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def recompute_adduct_networks(source, graph_store=None, accepted_only=True):
    """Filter network membership before scoring; retain input-row mappings for graphs.

    accepted_only=False explicitly enables exploratory accepted/excluded networks.
    Nonmembers retain missing scores, never their input's precomputed strength.
    """
    required = set(NETWORK_GROUPS + ["Metabolite name", "MS/MS spectrum", "Average Mz",
                                   "Overall Match Status"])
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Recomputing adduct networks requires: " + ", ".join(sorted(missing)))
    # Import the existing implementation, preserving its cosine/neutral-loss score
    # and >0.1 edge threshold, rather than reimplementing the network algorithm.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from PublicScripts.Lipids.LipidFunctions import network_analysis

    frame = source.copy().reset_index(drop=True)
    for column in NETWORK_GROUPS + ["Metabolite name", "MS/MS spectrum"]:
        frame[column] = frame[column].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.split().str.join(" ")
    frame["Dataset"] = frame["Dataset"].str.casefold().map({d.casefold(): d for d in datasets})
    frame["Polarity"] = frame["Polarity"].str.casefold()
    canonical = {name.casefold(): name for name in sorted(frame["Ontology"].unique(), reverse=True)}
    frame["Ontology"] = frame["Ontology"].str.casefold().map(canonical)
    frame["Average Mz"] = pd.to_numeric(frame["Average Mz"], errors="coerce")
    frame["Network Lipid ID"] = frame["Metabolite name"].str.split("|", n=1).str[0].str.strip().str.casefold()
    usable = (frame["Dataset"].notna() & frame["Ontology"].ne("")
              & frame["Polarity"].isin(POLARITY_COLORS) & frame["Adduct type"].ne("")
              & frame["Network Lipid ID"].ne("") & frame["MS/MS spectrum"].ne("")
              & np.isfinite(frame["Average Mz"]) & frame["Average Mz"].gt(0))
    accepted = frame["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    members = accepted if accepted_only else accepted | frame["Overall Match Status"].isin(
        [False, "FALSE", "False", "false", 0, "0"])
    if "Validation Status" in frame:
        members &= frame["Validation Status"].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.strip().str.casefold().ne("provisional")
    usable &= members
    frame["Class Strength"] = np.nan
    frame = frame.drop(columns=["Class Strength Percentile", "Max-Normalized Class Strength",
                                "Square Root Class Strength"], errors="ignore")
    frame["Network grouping"] = "Dataset / Ontology / Polarity / Adduct type"
    frame["Network population"] = "Accepted only" if accepted_only else "Accepted + excluded"
    frame["Network status"] = "Missing/invalid network inputs"
    frame.loc[~members, "Network status"] = "Not in network population"
    audit = []
    groups = list(frame.loc[usable].groupby(NETWORK_GROUPS, sort=True))
    for number, (key, group) in enumerate(groups, 1):
        row = dict(zip(NETWORK_GROUPS, key))
        row["Network population"] = "Accepted only" if accepted_only else "Accepted + excluded"
        tolerance = 0.05
        if "Tolerance" in group:
            values = pd.to_numeric(group["Tolerance"], errors="coerce")
            if not (np.isfinite(values) & values.gt(0)).all() or values.nunique() != 1:
                raise ValueError(f"{key}: expected one finite positive MS/MS Tolerance")
            tolerance = float(values.iloc[0])
        row.update({"Feature rows": len(group), "Unique lipids": group["Network Lipid ID"].nunique(),
                    "Tolerance": tolerance, "Edges": 0, "Cross-adduct edges": 0})
        if row["Unique lipids"] < 3:
            row["Status"] = "Skipped: fewer than 3 distinct lipids"
            frame.loc[group.index, "Network status"] = row["Status"]
            audit.append(row)
            continue
        print(f"Adduct network {number}/{len(groups)}: {key}, {len(group)} spectra", flush=True)
        scored, graph = network_analysis(group.reset_index(drop=True), mz_tol=tolerance,
                                         header="Class_", add_all_cols=True)
        if graph_store is not None:
            for node in graph.nodes():
                graph.nodes[node]["Input Row"] = group.index[int(node.split("~", 1)[0])]
            graph_store[key] = graph
        cross_edges = sum(graph.nodes[a]["Adduct type"] != graph.nodes[b]["Adduct type"]
                          for a, b in graph.edges)
        if cross_edges:
            raise ValueError(f"{key}: network contains cross-adduct edges")
        strength = scored["Class_strength"].to_numpy(dtype=float)
        if not (np.isfinite(strength) & (strength >= 0)).all():
            raise ValueError(f"{key}: invalid recomputed network strengths")
        frame.loc[group.index, "Class Strength"] = strength
        frame.loc[group.index, "Network status"] = "Recomputed within adduct"
        row.update({"Edges": graph.number_of_edges(), "Status": "Recomputed within adduct"})
        audit.append(row)
    return frame.drop(columns="Network Lipid ID"), pd.DataFrame(audit)


def normalization_audit(points):
    report = points.groupby(NETWORK_GROUPS, sort=True).agg(
        Lipids=("Lipid ID", "size"), Maximum=("Class Strength", "max"),
        **{"Normalized maximum": (NORMALIZED, "max")}).reset_index()
    report["Expected maximum"] = np.where(report["Maximum"].gt(0), 1.0, 0.0)
    report["Pass"] = np.isclose(report["Normalized maximum"], report["Expected maximum"])
    if not report["Pass"].all():
        raise ValueError("Within-adduct normalization maximum is not 1 (or 0 for an all-zero group)")
    return report


def comparable_points(points, min_lipids=MIN_LIPIDS, dataset_order=None, remove_zero_groups=True):
    """Require lipids in every dataset and the minimum in the largest dataset."""
    if min_lipids < 2:
        raise ValueError("min_lipids must be at least 2")
    required_datasets = datasets if dataset_order is None else dataset_order
    rows = []
    keep = []
    for (panel, polarity), group in points.groupby(["Plot Class", "Polarity"], sort=True):
        # Remove edgeless original-class/adduct groups before counting the panel.
        zero_counts = pd.Series(dtype=int)
        if remove_zero_groups:
            zero_counts = group.loc[group["Normalization maximum"].eq(0)].groupby("Dataset").size()
            group = group.loc[group["Normalization maximum"].gt(0)]
        counts = group.groupby("Dataset")["Lipid ID"].nunique()
        maximum_count = max((int(counts.get(d, 0)) for d in required_datasets), default=0)
        missing_datasets = [d for d in required_datasets if counts.get(d, 0) == 0]
        supported = maximum_count >= min_lipids and not missing_datasets
        reasons = []
        if maximum_count < min_lipids:
            reasons.append(f"Panel maximum across datasets is {maximum_count} eligible IDs; requires at least {min_lipids}")
        if missing_datasets:
            reasons.append("No eligible lipids in: " + ", ".join(missing_datasets))
        for dataset in required_datasets:
            retained = supported
            status = "Retained" if supported else "; ".join(reasons)
            rows.append({"Plot Class": panel, "Polarity": polarity, "Dataset": dataset,
                         "Lipids": counts.get(dataset, 0), "Minimum lipids": min_lipids,
                         "Maximum lipids across datasets": maximum_count,
                         "Zero-strength lipids removed": zero_counts.get(dataset, 0),
                         "Retained": retained, "Status": status})
            if retained:
                keep.extend(group.loc[group["Dataset"].eq(dataset)].index)
    return points.loc[keep].copy(), pd.DataFrame(rows)


def prepare_strength(source, all_ids=True, excluded_only=False, allow_empty=False):
    required = {"Dataset", "Ontology", "Polarity", "Metabolite name",
                "Overall Match Status", "Class Strength", "Adduct type"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))
    frame = source.copy()
    for column in ["Dataset", "Ontology", "Polarity", "Metabolite name", "Adduct type"]:
        frame[column] = frame[column].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.split().str.join(" ")
    frame["Dataset"] = frame["Dataset"].str.casefold().map({d.casefold(): d for d in datasets})
    frame["Polarity"] = frame["Polarity"].str.casefold()
    canonical = {name.casefold(): name for name in sorted(frame["Ontology"].unique(), reverse=True)}
    frame["Ontology"] = frame["Ontology"].str.casefold().map(canonical)
    frame["Lipid Name"] = frame["Metabolite name"].str.split("|", n=1).str[0].str.strip()
    frame["Lipid ID"] = frame["Ontology"].str.casefold() + "|" + frame["Lipid Name"].str.casefold()
    frame["Accepted"] = frame["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    failed = frame["Overall Match Status"].isin([False, "FALSE", "False", "false", 0, "0"])
    frame["Class Strength"] = pd.to_numeric(frame["Class Strength"], errors="coerce")
    valid = (frame["Dataset"].notna() & frame["Ontology"].ne("") & frame["Lipid Name"].ne("")
             & frame["Polarity"].isin(POLARITY_COLORS) & frame["Adduct type"].ne("")
             & (frame["Accepted"] | failed)
             & np.isfinite(frame["Class Strength"]) & frame["Class Strength"].ge(0))
    if "Validation Status" in frame:
        valid &= frame["Validation Status"].fillna("").astype(str).str.strip().str.casefold().ne("provisional")
    frame = frame.loc[valid].copy()
    frame["Plot Class"] = frame["Ontology"]
    classes = frame.loc[frame["Accepted"], "Plot Class"].unique()
    frame = frame.loc[frame["Plot Class"].isin(classes)].copy()
    frame["Class Group"] = frame["Plot Class"]
    frame["Plot Class"] = frame["Class Group"] + "\n" + frame["Adduct type"]
    accepted_only = not all_ids and not excluded_only
    if excluded_only:
        frame = frame.loc[~frame["Accepted"]].copy()
    elif accepted_only:
        frame = frame.loc[frame["Accepted"]].copy()
    keys = ["Plot Class", "Class Group", "Ontology", "Polarity", "Adduct type", "Lipid ID", "Dataset"]
    points = frame.groupby(keys, sort=True).agg(
        **{"Lipid Name": ("Lipid Name", "first"),
           "Class Strength": ("Class Strength", "median"),
           "Feature rows": ("Class Strength", "size"), "Accepted rows": ("Accepted", "sum")}
    ).reset_index()
    for _, point in points.loc[points["Feature rows"].gt(1)].iterrows():
        print(f"Merged {point['Feature rows']} rows for {point['Lipid Name']} in {point['Dataset']} / "
              f"{point['Polarity']} / {point['Adduct type']} using median Class Strength {point['Class Strength']}", flush=True)
    points["Normalization population"] = ("Excluded only" if excluded_only else
                                          "Accepted only" if accepted_only else "Accepted + excluded")
    points["Normalization unit"] = "Per-lipid-adduct median strength; maximum within adduct"
    points["Normalization maximum"] = points.groupby(
        NETWORK_GROUPS)["Class Strength"].transform("max")
    maximum = points["Normalization maximum"]
    points[NORMALIZED] = points["Class Strength"] / maximum.where(maximum.gt(0))
    points.loc[maximum.eq(0), NORMALIZED] = 0.0
    if points.empty and not allow_empty:
        raise ValueError("No eligible normalized class-strength scores")
    points["Excluded rows"] = points["Feature rows"] - points["Accepted rows"]
    normalization_audit(points)
    return points


def strength_summary(points, reference=None):
    rows = []
    layout = points if reference is None else reference
    for (lipid_class, polarity), _ in layout.groupby(["Plot Class", "Polarity"], sort=True):
        group = points.loc[points["Plot Class"].eq(lipid_class) & points["Polarity"].eq(polarity)]
        for dataset in datasets:
            values = group.loc[group["Dataset"].eq(dataset), NORMALIZED]
            rows.append({"Plot Class": lipid_class, "Polarity": polarity, "Dataset": dataset,
                         "Lipids": len(values), "Mean": values.mean(), "Median": values.median(),
                         "Q1": values.quantile(0.25), "Q3": values.quantile(0.75)})
    return pd.DataFrame(rows)


def compare_with_otot(summary):
    reference = summary.loc[summary["Dataset"].eq("OTOT")].drop(columns="Dataset").rename(
        columns={name: "OT/OT " + name for name in ["Lipids", "Mean", "Median", "Q1", "Q3"]})
    comparisons = summary.loc[summary["Dataset"].ne("OTOT")].merge(
        reference, on=["Plot Class", "Polarity"], how="left", validate="many_to_one")
    comparisons["Median difference (dataset - OT/OT)"] = comparisons["Median"] - comparisons["OT/OT Median"]
    comparisons["Mean difference (dataset - OT/OT)"] = comparisons["Mean"] - comparisons["OT/OT Mean"]
    comparisons["Status"] = np.where(
        comparisons["Lipids"].eq(0) | comparisons["OT/OT Lipids"].eq(0),
        "Unavailable: dataset or OT/OT has no scores",
        np.where(comparisons["Lipids"].lt(2) | comparisons["OT/OT Lipids"].lt(2),
                 "Descriptive only: singleton dataset or OT/OT", "Descriptive comparison"))
    return comparisons


def instrument_statistics(points, reference=None, all_pairs=False, value_column=NORMALIZED):
    """Lipid-clustered mean contrasts; Holm over all planned comparisons supplied."""
    rows = []
    pairs = list(combinations(datasets, 2)) if all_pairs else [("OTOT", d) for d in datasets[1:]]
    difference_column = "Mean difference (dataset - reference)" if all_pairs else "Mean difference (dataset - OT/OT)"
    layout = points if reference is None else reference
    for (lipid_class, polarity), _ in layout.groupby(["Plot Class", "Polarity"], sort=True):
        group = points.loc[points["Plot Class"].eq(lipid_class) & points["Polarity"].eq(polarity)]
        counts = group.groupby("Dataset").size()
        instruments = [d for d in datasets if counts.get(d, 0) >= 2]
        eligible = group.loc[group["Dataset"].isin(instruments)]
        means, covariance = None, None
        clusters = eligible["Lipid ID"].nunique()
        df = clusters - 1
        if len(instruments) >= 2 and clusters >= 2:
            design = np.column_stack([eligible["Dataset"].eq(d).to_numpy(dtype=float) for d in instruments])
            ids = pd.factorize(eligible["Lipid ID"])[0]
            result = OLS(eligible[value_column].to_numpy(dtype=float), design).fit()
            robust = result.get_robustcov_results(cov_type="cluster", groups=ids,
                                                use_correction=True, df_correction=True, use_t=True)
            means, covariance = np.asarray(robust.params), np.asarray(robust.cov_params())
        for baseline, dataset in pairs:
            row = {"Plot Class": lipid_class, "Polarity": polarity, "Dataset": dataset,
                   difference_column: np.nan, "Standard error": np.nan,
                   "CI95 low": np.nan, "CI95 high": np.nan, "p-value": np.nan,
                   "Lipid clusters": clusters, "DF": df, "Status": "OK"}
            if all_pairs:
                row["Reference Dataset"] = baseline
            if counts.get(dataset, 0) == 0 or counts.get(baseline, 0) == 0:
                row["Status"] = "Unavailable: dataset or reference has no scores"
            elif dataset not in instruments or baseline not in instruments:
                row["Status"] = "Unavailable: dataset or reference has fewer than two scores"
            elif means is None or not np.isfinite(covariance).all():
                row["Status"] = "Unavailable: covariance could not be estimated"
            else:
                i, j = instruments.index(dataset), instruments.index(baseline)
                difference = float(means[i] - means[j])
                variance = covariance[i, i] + covariance[j, j] - 2 * covariance[i, j]
                tolerance = 64 * np.finfo(float).eps * max(
                    abs(covariance[i, i]), abs(covariance[j, j]), np.finfo(float).tiny)
                constant = all(np.ptp(group.loc[group["Dataset"].eq(d), value_column]) == 0
                               for d in [dataset, baseline])
                row[difference_column] = difference
                if constant or not np.isfinite(variance) or variance <= tolerance:
                    row["Status"] = "Unavailable: no independent variation in the contrast"
                else:
                    se = float(np.sqrt(variance))
                    critical = t_distribution.ppf(0.975, df)
                    row.update({"Standard error": se, "CI95 low": difference - critical * se,
                                "CI95 high": difference + critical * se,
                                "p-value": 2 * t_distribution.sf(abs(difference / se), df)})
            rows.append(row)
    table = pd.DataFrame(rows)
    pvalues = table["p-value"].to_numpy()
    order = np.argsort(np.where(np.isfinite(pvalues), pvalues, 1.0))
    adjusted = np.empty(len(pvalues))
    adjusted[order] = np.minimum(1, np.maximum.accumulate(
        np.where(np.isfinite(pvalues[order]), pvalues[order], 1.0) * np.arange(len(pvalues), 0, -1)))
    adjusted[~np.isfinite(pvalues)] = np.nan
    table["Holm-adjusted p-value"] = adjusted
    table["Significant"] = np.isfinite(adjusted) & (adjusted < 0.05)
    return table


def panel_letter(index):
    """Return zero-based panel letters, continuing with AA after Z."""
    label = ""
    while index >= 0:
        index, remainder = divmod(index, 26)
        label = chr(ord("A") + remainder) + label
        index -= 1
    return label


def plot_strength(points, all_ids=True, statistics=None, excluded_only=False, reference=None,
                  value_column=NORMALIZED, ylabel=None, normalized_axis=True):
    layout = points if reference is None else reference
    nonzero = points.loc[points[value_column].gt(0)] if normalized_axis else points
    classes = {mode: set(layout.loc[layout["Polarity"].eq(mode), "Plot Class"])
               & set(nonzero.loc[nonzero["Polarity"].eq(mode), "Plot Class"])
               for mode in POLARITY_COLORS}
    common = sorted(classes["positive"] & classes["negative"])
    slots = {mode: common + sorted(names - set(common)) for mode, names in classes.items()}
    figures = {}
    pages = []
    for mode, color in POLARITY_COLORS.items():
        count = max(1, int(np.ceil(len(slots[mode]) / MAX_PANELS_PER_PAGE)))
        per_page = MAX_PANELS_PER_PAGE
        for page in range(count):
            pages.append((mode, color, page, count,
                          slots[mode][page * per_page:(page + 1) * per_page]))
    for mode, color, page, page_count, mode_slots in pages:
        panel_offset = page * MAX_PANELS_PER_PAGE
        significant = None
        max_brackets = 0
        if statistics is not None:
            significant = statistics.loc[statistics["Polarity"].eq(mode)
                                         & statistics["Plot Class"].isin(mode_slots)
                                         & statistics["Significant"]].copy()
            if not significant.empty:
                max_brackets = int(significant.groupby("Plot Class").size().max())
        scale = 1.0
        if not normalized_axis:
            values = points.loc[points["Polarity"].eq(mode) & points["Plot Class"].isin(mode_slots), value_column]
            scale = max(1.0, float(values.max())) if len(values) else 1.0
        upper = scale * (1.18 if not max_brackets else 1.58 + 0.16 * (max_brackets - 1))
        columns = 4
        rows = max(1, int(np.ceil(len(mode_slots) / columns)))
        figure_size = (FIGURE_SIZE[0], FIGURE_SIZE[1] * rows / 5)
        fig, axes = plt.subplots(rows, columns, figsize=figure_size,
                                 squeeze=False, layout="constrained")
        rng = np.random.default_rng(42)
        for index, ax in enumerate(axes.flat):
            if index >= len(mode_slots):
                ax.set_axis_off()
                continue
            lipid_class = mode_slots[index]
            ax.text(0.03, 0.95, panel_letter(panel_offset + index), transform=ax.transAxes,
                    ha="left", va="top", fontsize=8, fontweight="bold")
            selected = points.loc[points["Polarity"].eq(mode) & points["Plot Class"].eq(lipid_class)]
            for position, dataset in enumerate(datasets):
                values = selected.loc[selected["Dataset"].eq(dataset), value_column].to_numpy()
                if len(values):
                    ax.scatter(position + rng.uniform(-0.12, 0.12, len(values)), values,
                               color=color, alpha=0.45, s=8, linewidths=0, rasterized=True)
                    ax.boxplot([values], positions=[position], widths=0.4, manage_ticks=False,
                               showfliers=False, patch_artist=True,
                               boxprops={"facecolor": to_rgba(color, 0.15), "edgecolor": color},
                               medianprops={"color": color, "linewidth": 1.3},
                               whiskerprops={"color": color}, capprops={"color": color})
                ax.text(position, 1.10 * scale, str(len(values)), ha="center", fontsize=8)
            if significant is not None:
                comparisons = significant.loc[significant["Plot Class"].eq(lipid_class)].copy()
                # Short brackets go below long brackets to keep their endpoints clear.
                comparisons["Span"] = [datasets.index(row["Dataset"])
                                       - datasets.index(row.get("Reference Dataset", "OTOT"))
                                       for _, row in comparisons.iterrows()]
                comparisons = comparisons.sort_values("Span", kind="stable")
                for level, (_, comparison) in enumerate(comparisons.iterrows()):
                    left = datasets.index(comparison.get("Reference Dataset", "OTOT"))
                    right = datasets.index(comparison["Dataset"])
                    height = scale * (1.38 + 0.16 * level)
                    ax.plot([left, left, right, right],
                            [height - 0.035 * scale, height, height, height - 0.035 * scale],
                            color="black", linewidth=0.7)
                    p = comparison["Holm-adjusted p-value"]
                    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*"
                    ax.text((left + right) / 2, height, stars,
                            ha="center", va="center", fontsize=8,
                            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.3})
            ax.set_title(lipid_class, fontsize=9)
            ax.set_ylim(-0.03 * scale, upper)
            if normalized_axis:
                ax.set_yticks([0, 0.5, 1])
            else:
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.set_xticks(range(len(datasets)), labels=[DATASET_LABELS[d] for d in datasets],
                          rotation=40, ha="right")
            ax.tick_params(axis="x", labelsize=7,
                           labelbottom=index + columns >= len(mode_slots))
            ax.tick_params(axis="y", labelsize=8)
            ax.set_xlim(-0.6, len(datasets) - 0.4)
            ax.spines[["top", "right"]].set_visible(False)
            if index % columns:
                ax.tick_params(labelleft=False)
        if not mode_slots:
            axes[0, 0].set_axis_on()
            axes[0, 0].text(0.5, 0.5, "No eligible scores", ha="center", transform=axes[0, 0].transAxes)
        population = "Excluded only" if excluded_only else "Accepted + excluded" if all_ids else "Accepted only"
        page_label = f" | Page {page + 1}/{page_count}" if page_count > 1 else ""
        fig.suptitle(f"{mode.title()} | {population}{page_label}", fontsize=10)
        fig.supylabel(ylabel or "Class Strength / Within-Dataset Class-Adduct Maximum", fontsize=10)
        if statistics is not None:
            fig.supxlabel("Holm-adjusted p: * <0.05, ** <0.01, *** <0.001", fontsize=8)
        key = f"{mode}_Page_{page + 1}" if page_count > 1 else mode
        figures[key] = fig
    return figures


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-lipids", type=int, default=MIN_LIPIDS,
                        help="Minimum for the maximum unique eligible lipid count across datasets; every dataset must have at least one (default: 10)")
    population = parser.add_mutually_exclusive_group()
    population.add_argument("--accepted-only", action="store_true")
    population.add_argument("--excluded-only", action="store_true",
                            help="Plot excluded IDs scored in exploratory accepted + excluded networks")
    population.add_argument("--all-ids", action="store_true",
                            help="Build and plot exploratory accepted + excluded networks")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    if args.min_lipids < 2:
        parser.error("--min-lipids must be at least 2")

    source = pd.read_csv(args.input) if args.input.suffix.lower() == ".csv" else pd.read_excel(args.input, sheet_name="Sheet1")

    source, _ = recompute_adduct_networks(source, accepted_only=not (args.all_ids or args.excluded_only))
    output_root = args.output_dir or args.input.parent / "Final Figures"
    output_root.mkdir(parents=True, exist_ok=True)
    reference = prepare_strength(source, allow_empty=True)
    populations = [("Accepted_IDs", False, False)]
    if args.excluded_only:
        populations = [("Excluded_IDs", True, True)]
    elif args.all_ids:
        populations = [("All_IDs", True, False)]
    for name, all_ids, excluded_only in populations:
        output_dir = output_root
        export_population(source, output_dir, args.input, all_ids, args.no_show, excluded_only, reference,
                          min_lipids=args.min_lipids)
    if not args.no_show:
        plt.show()


def export_population(source, output_dir, input_file, all_ids, no_show, excluded_only=False, reference=None,
                      min_lipids=MIN_LIPIDS):
    points = prepare_strength(source, all_ids=all_ids, excluded_only=excluded_only, allow_empty=True)
    points, _ = comparable_points(points, min_lipids=min_lipids)
    # Statistics and multiplicity correction use exactly the retained panels.
    reference = points
    output_dir.mkdir(parents=True, exist_ok=True)
    normalization_audit(points)
    if points.empty:
        print(f"No comparable panels in {output_dir}.")
        return
    summary = strength_summary(points, reference)
    comparisons = compare_with_otot(summary)
    pairwise = instrument_statistics(points, reference, all_pairs=True)
    figures = plot_strength(points, all_ids, pairwise, excluded_only, reference)
    for mode, figure in figures.items():
        polarity, separator, page = mode.partition("_Page_")
        basename = OUTPUT_BASENAMES[polarity]
        if excluded_only:
            basename += "_Excluded_IDs"
        elif all_ids:
            basename += "_All_IDs"
        if separator:
            basename += f"_Page_{page}"
        # Resolve the dense constrained layout before the first high-DPI export.
        figure.canvas.draw()
        for extension in ["png", "pdf"]:
            figure.savefig(output_dir / f"{basename}.{extension}", dpi=600)
        if no_show:
            plt.close(figure)
    print(comparisons[["Plot Class", "Polarity", "Dataset", "Lipids", "Median", "OT/OT Median",
                       "Median difference (dataset - OT/OT)", "Status"]].to_string(index=False, na_rep="N/A"))
    print(f"Saved class-strength comparison in {output_dir}")


if __name__ == "__main__":
    main()
