"""Plot SI_10 tail composition and SI_11 tail-product rarity versus dataset coverage.

Reuse accepted-ID selection, resolved chains, multiplicities, and deduplication
from Figure5. Save only PNG/PDF figures: SI_10 and SI_11 in Final Figures.

Plot lipid-level tail-product rarity against dataset coverage, stratified
by tail count. Use all tail frequencies (no display cutoff), averaged equally
across nonempty datasets, expressed as fractions. For each unique
class/chain composition with >=1 fatty-acyl tails, P = product(f_tail), with
repeated chains included repeatedly; rarity = -log10(P). Exclude the whole
lipid if any chain is LCB or O-/P-, rather than scoring a partial composition.
Frequencies retain the full annotated-tail denominator. Calculate component
frequencies and scores in memory. This is an exploratory
composition score, not a joint probability, abundance, or identification
confidence: chains are dependent and frequencies come from these same IDs.
Panels separate tail counts because products decrease as chain count increases.
Single-tail scores are -log10(f_tail). Violins show distributions with black
median marks; dots are unique lipids with visual-only jitter (seed 42).
Singleton/constant groups use a horizontal mark instead of a density estimate.
Match the normalized-intensity analysis: categorical coverage OLS mean tests,
CR1 covariance clustered by unique class/chain composition, G-1 denominator df,
and two-sided pairwise t contrasts. Each lipid appears once, so clusters are
singletons (equivalent to HC1); dataset fixed effects are not applicable to
these pooled composition scores. Holm-adjust overall tests across tail-count
panels and pairwise tests across all panel/coverage pairs (4 and 40 planned
tests with the default 1-4 tail panels and five datasets). Additional observed
tail counts extend both families. Missing/singleton strata are unavailable;
unavailable planned tests count as p=1 in correction. Stars show adjusted
p < .05/.01/.001/.0001; only significant comparisons receive brackets.
These are exploratory, conditional tests: shared tails and frequencies estimated
from these same IDs create between-lipid dependencies that lipid clustering
does not address. They do not establish independent instrument replication.
"""

import argparse
from itertools import combinations
from pathlib import Path
import sys

import matplotlib
if __name__ == "__main__" and "--no-show" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests

from Figure5 import (
    DATASETS, DATASET_LABELS, TOPDIR, frequency_summary, prepare_lipids,
    tail_composition_table,
)


INPUT_FILE = TOPDIR / "Combined_Fragment_Assignment_Colored.xlsx"
PRODUCT_SCORE = "Tail-product rarity (-log10)"
PRODUCT_COLORS = ["#D9EAF4", "#B8D8E8", "#8FC1DA", "#5FA6C8", "#347FA8"]
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})



def lipid_product_tables(lipids, occurrences, datasets):
    """Score pooled unique lipids using unfiltered mean tail-copy fractions."""
    lipids = lipids.loc[lipids["Dataset"].isin(datasets)]
    occurrences = occurrences.loc[occurrences["Dataset"].isin(datasets)]
    summary = frequency_summary(lipids, occurrences, ["Dataset"])
    fractions = summary.pivot(index="Tail", columns="Dataset", values="Tail composition (%)")
    fractions = fractions.reindex(columns=datasets) / 100
    frequencies = fractions.mean(axis=1).rename("Mean tail fraction").to_frame()
    frequencies["Nonempty datasets"] = fractions.notna().sum(axis=1)
    for dataset in datasets:
        frequencies[f"{dataset} tail fraction"] = fractions[dataset]
    pooled = lipids.drop_duplicates(["Ontology", "Species"])[["Ontology", "Species", "Chains"]].copy()
    coverage = lipids.groupby(["Ontology", "Species"])["Dataset"].nunique()
    pooled["Dataset Coverage"] = [coverage.loc[(row.Ontology, row.Species)] for row in pooled.itertuples()]
    pooled["Tail count"] = pooled["Chains"].map(len)
    components, products, scores, reasons = [], [], [], []
    for chains in pooled["Chains"]:
        values = frequencies["Mean tail fraction"].reindex(chains).to_numpy(dtype=float)
        reason = ""
        if not chains:
            reason = "No resolved tails"
        elif any(tail.startswith(("LCB ", "O-", "P-")) for tail in chains):
            reason = "Contains LCB or ether tail"
        elif not np.isfinite(values).all() or np.any(values <= 0) or np.any(values > 1):
            reason = "Unavailable or invalid tail frequency"
        components.append("; ".join(f"{value:.12g}" for value in values))
        products.append(np.prod(values) if not reason else np.nan)
        scores.append(-np.log10(values).sum() if not reason else np.nan)
        reasons.append(reason)
    pooled["Mean tail fractions (chain order)"] = components
    pooled["Tail frequency product"] = products
    pooled["Tail-product rarity (-log10)"] = scores
    pooled["Score exclusion reason"] = reasons
    pooled["Chains"] = pooled["Chains"].map("; ".join)
    return pooled.reset_index(drop=True), frequencies.rename_axis("Tail").reset_index()


def product_significance(pvalue):
    if not np.isfinite(pvalue) or pvalue >= .05:
        return ""
    return "****" if pvalue < .0001 else "***" if pvalue < .001 else "**" if pvalue < .01 else "*"


def lipid_product_statistics(scored, datasets):
    """Match intensity mean tests, with one observation per unique lipid."""
    eligible = scored.loc[scored["Score exclusion reason"].eq("")]
    if eligible.duplicated(["Ontology", "Species"]).any():
        raise ValueError("Product statistics require one observation per unique lipid")
    tail_counts = sorted({1, 2, 3, 4} | set(eligible["Tail count"]))
    overall, pairwise, summary = [], [], []
    for count in tail_counts:
        group = eligible.loc[eligible["Tail count"].eq(count)]
        sizes = group.groupby("Dataset Coverage").size()
        levels = sorted(sizes.index)
        n = len(group)
        result = {"Tail count": count, "Lipids": n,
                  "F statistic": np.nan, "Numerator df": np.nan, "Denominator df": n - 1,
                  "p-value": np.nan, "Status": "Unavailable: insufficient data or singular covariance"}
        fit = None
        columns = levels[1:]
        if len(levels) > 1 and all(sizes[level] >= 2 for level in levels) and group[PRODUCT_SCORE].nunique() > 1:
            design = np.column_stack([np.ones(n)] + [group["Dataset Coverage"].eq(level) for level in columns])
            fit = OLS(group[PRODUCT_SCORE].to_numpy(), design).fit().get_robustcov_results(
                cov_type="cluster", groups=np.arange(n), use_correction=True,
                df_correction=True, use_t=True)
            contrasts = np.eye(design.shape[1])[1:]
            covariance = contrasts @ fit.cov_params() @ contrasts.T
            if np.isfinite(covariance).all() and np.linalg.matrix_rank(covariance) == len(columns):
                test = fit.f_test(contrasts)
                result.update({"F statistic": float(np.asarray(test.fvalue).item()),
                               "Numerator df": test.df_num, "p-value": float(test.pvalue),
                               "Status": "Exploratory clustered mean test"})
        overall.append(result)
        for coverage in range(1, len(datasets) + 1):
            values = group.loc[group["Dataset Coverage"].eq(coverage), PRODUCT_SCORE]
            summary.append({"Tail count": count, "Dataset Coverage": coverage, "Lipids": len(values),
                            "Mean rarity": values.mean(), "Median rarity": values.median(),
                            "Q1": values.quantile(.25), "Q3": values.quantile(.75)})
        for left, right in combinations(range(1, len(datasets) + 1), 2):
            comparison = {"Tail count": count, "Coverage A": left, "Coverage B": right,
                          "Lipids A": int(sizes.get(left, 0)), "Lipids B": int(sizes.get(right, 0)),
                          "Mean difference (B - A)": np.nan, "p-value": np.nan,
                          "Status": "Unavailable: insufficient data or singular covariance"}
            if fit is not None and left in levels and right in levels:
                contrast = np.zeros(len(fit.params))
                for level, sign in [(left, -1), (right, 1)]:
                    if level in columns:
                        contrast[1 + columns.index(level)] = sign
                variance = float(contrast @ fit.cov_params() @ contrast)
                if np.isfinite(variance) and variance > 0:
                    test = fit.t_test(contrast)
                    comparison.update({"Mean difference (B - A)": float(np.asarray(test.effect).item()),
                                       "p-value": float(np.asarray(test.pvalue).item()),
                                       "Status": "Exploratory clustered mean contrast"})
            pairwise.append(comparison)
    overall, pairwise = pd.DataFrame(overall), pd.DataFrame(pairwise)
    for table in [overall, pairwise]:
        available = np.isfinite(table["p-value"])
        adjusted = multipletests(table["p-value"].fillna(1), method="holm")[1]
        table["Holm-adjusted p-value"] = np.where(available, adjusted, np.nan)
    pairwise["Significance"] = pairwise["Holm-adjusted p-value"].map(product_significance)
    return overall, pairwise, pd.DataFrame(summary)


def plot_lipid_products(scored, datasets, statistics, pairwise):
    """Violin panels by tail count, matching the normalized-intensity format."""
    eligible = scored.loc[scored["Score exclusion reason"].eq("")]
    tail_counts = statistics["Tail count"].tolist()
    significant = pairwise.loc[pairwise["Significance"].ne("")]
    bracket_rows = max(1, int(significant.groupby("Tail count").size().max())) if len(significant) else 1
    rows = (len(tail_counts) + 1) // 2
    fig = plt.figure(figsize=(6.0, 6.25 * rows / 2), layout="constrained")
    grid = fig.add_gridspec(rows, 2)
    axes = []
    rng = np.random.default_rng(42)
    maximum = max(1, eligible[PRODUCT_SCORE].max()) if len(eligible) else 1
    colors = PRODUCT_COLORS if len(datasets) <= 5 else plt.get_cmap("Blues")(np.linspace(.2, .75, len(datasets)))
    for index, count in enumerate(tail_counts):
        cell = grid[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=[.18 + .14 * bracket_rows, 1.55], hspace=.02)
        bracket_ax = fig.add_subplot(cell[0])
        ax = fig.add_subplot(cell[1], sharex=axes[0] if axes else None, sharey=axes[0] if axes else None)
        axes.append(ax)
        ax.tick_params(labelbottom=index // 2 == rows - 1, labelleft=index % 2 == 0)
        group = eligible.loc[eligible["Tail count"].eq(count)]
        for coverage, color in zip(range(1, len(datasets) + 1), colors):
            values = group.loc[group["Dataset Coverage"].eq(coverage), PRODUCT_SCORE].to_numpy()
            if len(values) > 1 and np.ptp(values) > 0:
                violin = ax.violinplot([values], positions=[coverage], widths=.7,
                                      showmedians=True, showextrema=False)
                for body in violin["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(color)
                    body.set_alpha(.6)
                violin["cmedians"].set_color("black")
            elif len(values):
                ax.hlines(values[0], coverage - .2, coverage + .2, color="black", linewidth=1)
            ax.scatter(coverage + rng.uniform(-.12, .12, len(values)), values,
                       s=5, color="#347FA8", alpha=.35, linewidths=0, zorder=3)
            ax.text(coverage, .96, f"n={len(values)}", ha="center", va="center", fontsize=7,
                    transform=ax.get_xaxis_transform())
        result = statistics.loc[statistics["Tail count"].eq(count)].iloc[0]
        p = result["Holm-adjusted p-value"]
        p_text = f"p\u2090={p:.2g}" if np.isfinite(p) else "p\u2090=N/A"
        ax.set_title(f"{chr(65 + index)}  {count} Tail{'s' if count != 1 else ''}",
                             loc="left", fontweight="bold")
        ax.set_title(p_text, loc="right", fontsize=7)
        pairs = significant.loc[significant["Tail count"].eq(count)].copy()
        pairs["Span"] = pairs["Coverage B"] - pairs["Coverage A"]
        for row, (_, pair) in enumerate(pairs.sort_values(["Span", "Coverage A"]).iterrows()):
            left, right = pair["Coverage A"], pair["Coverage B"]
            height = row + .25
            bracket_ax.plot([left, left, right, right], [height, height + .15, height + .15, height],
                            color="black", linewidth=.7)
            bracket_ax.text((left + right) / 2, height + .16, pair["Significance"],
                            ha="center", va="bottom", fontsize=8)
        if pairs.empty:
            available = pairwise.loc[pairwise["Tail count"].eq(count), "Holm-adjusted p-value"].notna().any()
            bracket_ax.text(.5, .5, "No Significant Pairs" if available else "Pairwise Tests Unavailable",
                            transform=bracket_ax.transAxes, ha="center", va="center", fontsize=7)
        bracket_ax.set_xlim(.5, len(datasets) + .5)
        bracket_ax.set_ylim(0, bracket_rows + .3)
        bracket_ax.set_axis_off()
        ax.set_xticks(range(1, len(datasets) + 1))
        ax.set_xlim(0.5, len(datasets) + 0.5)
        ax.set_ylim(-maximum * .03, maximum * 1.15)
        ax.spines[["top", "right"]].set_visible(False)
    fig.supxlabel("Datasets identifying the lipid composition", fontsize=9)
    fig.supylabel("Tail-product rarity " + r"$-\log_{10}(P)$", fontsize=9)
    return fig


def plot_tail_overview(percentages, population="Accepted IDs"):
    """Show unmodified tail percentages on one shared, zero-based color scale."""
    rows, columns = percentages.shape
    fig, ax = plt.subplots(figsize=(6.0, 6.25), layout="constrained")
    if percentages.empty:
        ax.text(0.5, 0.5, "No tails meet the display threshold",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        values = percentages.to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        norm = Normalize(vmin=0, vmax=max(1.0, float(finite.max()) if finite.size else 0))
        cmap = plt.get_cmap("Blues").copy()
        cmap.set_bad("#E5E5E5")
        image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(columns),
                      labels=[DATASET_LABELS[d] for d in percentages.columns])
        ax.set_yticks(range(rows), labels=percentages.index)
        ax.tick_params(length=0, labelsize=8, pad=5)
        ax.set_xlabel("Acquisition Method", fontsize=9, labelpad=7)
        ax.set_ylabel("Tail Identity", fontsize=9, labelpad=7)
        ax.set_xticks(np.arange(columns + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(rows + 1) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="minor", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for row in range(rows):
            for column in range(columns):
                value = values[row, column]
                label = f"{value:.1f}%" if np.isfinite(value) else "N/A"
                rgba = cmap(norm(value)) if np.isfinite(value) else (0.9, 0.9, 0.9, 1)
                luminance = np.dot(rgba[:3], [0.299, 0.587, 0.114])
                ax.text(column, row, label, ha="center", va="center",
                        fontsize=8, color="white" if luminance < 0.5 else "#222222")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.04, shrink=0.7)
        colorbar.set_label("Annotated Tail Copies (%)", fontsize=9, labelpad=7)
        colorbar.ax.tick_params(labelsize=8, width=0.7)
        colorbar.outline.set_linewidth(0.7)
        colorbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_title(f"Tail Composition ({population})", fontsize=10, pad=10)
    return fig



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--datasets", nargs="+", choices=DATASETS, default=DATASETS)
    parser.add_argument("--min-tail-percent", type=float, default=1.0)
    parser.add_argument("--all-ids", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    if not np.isfinite(args.min_tail_percent) or not 0 <= args.min_tail_percent <= 100:
        parser.error("--min-tail-percent must be between 0 and 100")
    datasets = list(dict.fromkeys(args.datasets))
    if len(datasets) < 2:
        parser.error("Select at least two datasets for coverage comparison")
    lipids, occurrences, _ = prepare_lipids(args.input, datasets, args.all_ids)
    summary = frequency_summary(lipids, occurrences, ["Dataset"])
    percentages = tail_composition_table(summary, datasets, args.min_tail_percent)
    output = args.output_dir or args.input.parent / "Final Figures"
    output.mkdir(parents=True, exist_ok=True)
    overview = plot_tail_overview(percentages, "Accepted + Unvalidated IDs" if args.all_ids else "Accepted IDs")
    for extension in ["png", "pdf"]:
        overview.savefig(output / f"SI_10.{extension}", dpi=600)
    scored, frequencies = lipid_product_tables(lipids, occurrences, datasets)
    product_statistics, product_pairwise, product_summary = lipid_product_statistics(scored, datasets)
    product_fig = plot_lipid_products(scored, datasets, product_statistics, product_pairwise)
    for extension in ["png", "pdf"]:
        product_fig.savefig(output / f"SI_11.{extension}", dpi=600)
    print("Tail-product score accounting:")
    print(scored["Score exclusion reason"].replace("", "Scored lipid composition").value_counts().to_string())
    print(product_statistics.to_string(index=False))
    print("Tail-product mean tests are exploratory; stars are Holm-adjusted. Black marks are medians.")
    print(f"Saved comparison to {output}")
    if args.no_show:
        plt.close(overview)
        plt.close(product_fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
