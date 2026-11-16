"""Compare normalized class-network strength with confirmed dataset coverage.

Build adduct-specific networks from accepted, non-provisional spectra only and normalize lipid
strengths using SI_7_8. Coverage counts distinct manuscript datasets confirming
each original-ontology/sum-composition identity, including confirmations without
a usable network score. Score each identity by first averaging its normalized
adduct/polarity scores within each dataset, then averaging across scored datasets.
Each lipid contributes once. Missing scores are omitted, never replaced by zero;
valid all-zero networks retain zero, as in SI_7_8. No panel-size filter is applied.
Violins show score distributions and medians; jittered dots show every lipid
in the same muted blue. Black dots and capped bars show medians and pointwise
95% percentile-bootstrap confidence intervals (9,999 resamples, seed 42).
Resample lipid identities within each dataset/coverage group, using one score
per identity; the first figure's Combined panel resamples the already averaged
lipid scores. The second figure's Combined panel pools coverage within each
analyzer before resampling its lipid scores. Intervals
condition on the observed networks, normalization, and coverage assignments;
they do not quantify network-estimation or instrument-run uncertainty. Empty
and singleton groups have no interval; constant groups have degenerate bounds.
Save SI_5 and SI_6 as PNG/PDF figures only.
The first six-panel figure puts Combined in panel A, followed by the five
datasets in panels B-F. It compares coverage-group means within each dataset
and in Combined using robust ANOVA-style F tests and two-sided t-tests of mean contrasts.
Each lipid contributes once per panel; singleton lipid clusters give HC1-equivalent
covariance with n-1 denominator degrees of freedom. Combined uses the already
averaged lipid scores, not pooled lipid-dataset observations. Holm-adjust six
overall tests together and 60 planned pairwise tests together, separately from
the second figure. Retain shared axes and a 6.5 x 6.75-inch canvas.
A second 6.5 x 6.75-inch figure compares analyzers in panel order Combined,
5, 4, 3, 2, 1 in lettered panels A-F, with shared axes. Panel A (Combined) pools
all coverage levels, retaining one within-dataset score per lipid per analyzer.
Each lipid contributes its within-dataset mean score once per analyzer.
Use categorical-analyzer OLS mean comparisons with CR1 covariance clustered
by lipid identity and G-1 denominator degrees of freedom. The overall robust
ANOVA-style F test compares analyzer means; pairwise tests are two-sided t
contrasts from the same model. Shared lipid identities stay in the same cluster.
These overall F tests jointly compare means, rather than testing equal variances;
the pairwise t-tests compare individual mean differences.
Combined compares marginal analyzer means across all coverage levels, without
coverage adjustment. Tests concern means; plotted black marks remain medians.
Holm-adjust the six overall tests together and all 60 planned pairwise tests
together, counting unavailable tests as p=1 in each correction family.
Groups with fewer than two scored lipids are descriptive only. P-values are
exploratory: shared networks induce dependence and these datasets do not provide
independent instrument-run replication. Save only the comparison PNG/PDF figures.
"""

import argparse
from itertools import combinations
from pathlib import Path
import sys

import matplotlib
if __name__ == "__main__" and "--no-show" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests

from Figure3 import DATASET_LABELS, coverage_colors, datasets, topdir
from SI_8_9 import NORMALIZED, prepare_strength, recompute_adduct_networks


INPUT_FILE = Path(topdir) / "Combined_Fragment_Assignment_Colored.xlsx"
SCORE = "Mean normalized score"
BOOTSTRAPS = 9999

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def prepare_comparison(source, points):
    """Join accepted-ID coverage to one equal-dataset-weight score per lipid."""
    required = {"Dataset", "Ontology", "Metabolite name", "Overall Match Status"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))
    frame = source.copy()
    for column in ["Dataset", "Ontology", "Metabolite name"]:
        frame[column] = frame[column].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.split().str.join(" ")
    frame["Dataset"] = frame["Dataset"].str.casefold().map({d.casefold(): d for d in datasets})
    frame["Lipid Name"] = frame["Metabolite name"].str.split("|", n=1).str[0].str.strip()
    frame["Lipid ID"] = frame["Ontology"].str.casefold() + "|" + frame["Lipid Name"].str.casefold()
    valid = (frame["Dataset"].notna() & frame["Ontology"].ne("") & frame["Lipid Name"].ne("")
             & frame["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"]))
    if "Validation Status" in frame:
        valid &= frame["Validation Status"].fillna("").astype(str).str.strip().str.casefold().ne("provisional")
    coverage = frame.loc[valid].groupby("Lipid ID", sort=True).agg(
        **{"Lipid Name": ("Lipid Name", "first"), "Ontology": ("Ontology", "first"),
           "Dataset coverage": ("Dataset", "nunique")})
    # Coverage is independent of score availability. Give each scored dataset
    # equal weight even when datasets have different numbers of adducts.
    scored = points.loc[np.isfinite(points[NORMALIZED])].copy()
    per_dataset = scored.groupby(["Lipid ID", "Dataset"])[NORMALIZED].mean().reset_index()
    scores = per_dataset.groupby("Lipid ID").agg(
        **{SCORE: (NORMALIZED, "mean"), "Scored datasets": ("Dataset", "nunique")})
    paired = coverage.join(scores, how="left").reset_index()
    if not paired[SCORE].notna().any():
        raise ValueError("No confirmed lipids have eligible normalized scores")
    return paired


def prepare_dataset_comparison(paired, points):
    """Keep dataset-specific scores, joining coverage across all five datasets."""
    scored = points.loc[np.isfinite(points[NORMALIZED])]
    per_dataset = scored.groupby(["Lipid ID", "Dataset"])[NORMALIZED].mean().rename(SCORE).reset_index()
    return per_dataset.merge(
        paired[["Lipid ID", "Lipid Name", "Ontology", "Dataset coverage"]],
        on="Lipid ID", how="inner", validate="many_to_one")


def compare_coverage_scores(by_dataset, coverage):
    """Robust analyzer mean tests within a coverage stratum, or pooled for None."""
    if coverage is not None and coverage not in range(1, len(datasets) + 1):
        raise ValueError("Coverage must be between 1 and 5")
    selected = by_dataset.loc[np.isfinite(by_dataset[SCORE])].copy()
    if coverage is not None:
        selected = selected.loc[selected["Dataset coverage"].eq(coverage)].copy()
    if selected.duplicated(["Lipid ID", "Dataset"]).any():
        raise ValueError("Each lipid must contribute exactly once per dataset")
    counts = selected.groupby("Dataset").size()
    levels = [name for name in datasets if counts.get(name, 0) >= 2]
    # Missing and singleton analyzers remain descriptive, without preventing
    # comparisons between other analyzers that have estimable variability.
    modeled = selected.loc[selected["Dataset"].isin(levels)]
    clusters = modeled["Lipid ID"].nunique()
    label = coverage if coverage is not None else "All"
    overall = {"Dataset coverage": label, "Observations": len(modeled),
               "Lipid clusters": clusters, "Analyzers tested": "; ".join(levels),
               "F statistic": np.nan, "Numerator df": np.nan,
               "Denominator df": clusters - 1, "p-value": np.nan,
               "Status": "Unavailable: insufficient data or singular covariance"}
    fit = None
    if len(levels) > 1 and clusters > 1 and modeled[SCORE].nunique() > 1:
        design = np.column_stack([np.ones(len(modeled))] + [
            modeled["Dataset"].eq(name).to_numpy(dtype=float) for name in levels[1:]])
        if len(modeled) > design.shape[1]:
            fit = OLS(modeled[SCORE].to_numpy(), design).fit().get_robustcov_results(
                cov_type="cluster", groups=modeled["Lipid ID"].to_numpy(),
                use_correction=True, df_correction=True, use_t=True)
            contrasts = np.eye(len(levels))[1:]
            covariance = contrasts @ fit.cov_params() @ contrasts.T
            if np.isfinite(covariance).all() and np.linalg.matrix_rank(covariance) == len(levels) - 1:
                test = fit.f_test(contrasts)
                overall.update({"F statistic": float(np.asarray(test.fvalue).item()),
                                "Numerator df": test.df_num, "p-value": float(test.pvalue),
                                "Status": "Exploratory clustered analyzer mean test"})
    rows = []
    for left, right in combinations(datasets, 2):
        x = selected.loc[selected["Dataset"].eq(left)].set_index("Lipid ID")[SCORE]
        y = selected.loc[selected["Dataset"].eq(right)].set_index("Lipid ID")[SCORE]
        row = {"Dataset coverage": label, "Shared lipids": len(x.index.intersection(y.index)),
               "Reference Dataset": left, "Dataset": right,
               "Reference lipids": len(x), "Lipids": len(y),
               "Reference mean": x.mean(), "Mean": y.mean(),
               "Reference median": x.median(), "Median": y.median(),
               "Mean difference (Dataset - Reference)": np.nan,
               "Standard error": np.nan, "t statistic": np.nan,
               "Denominator df": clusters - 1, "p-value": np.nan,
               "Status": "Unavailable: insufficient data or singular covariance"}
        if fit is not None and left in levels and right in levels:
            contrast = np.zeros(len(levels))
            for name, sign in [(left, -1), (right, 1)]:
                if name != levels[0]:
                    contrast[levels.index(name)] = sign
            variance = float(contrast @ fit.cov_params() @ contrast)
            if np.isfinite(variance) and variance > 0:
                test = fit.t_test(contrast)
                row.update({"Mean difference (Dataset - Reference)": float(np.asarray(test.effect).item()),
                            "Standard error": float(np.asarray(test.sd).item()),
                            "t statistic": float(np.asarray(test.tvalue).item()),
                            "p-value": float(np.asarray(test.pvalue).item()),
                            "Status": "Exploratory clustered analyzer mean contrast"})
        rows.append(row)
    return selected, pd.DataFrame(rows), overall


def analyzer_comparisons(by_dataset):
    """Correct all planned panel tests together, including unavailable tests."""
    comparisons = []
    overall_rows = []
    for coverage in [None, *range(len(datasets), 0, -1)]:
        selected, pairwise, overall = compare_coverage_scores(by_dataset, coverage)
        comparisons.append((selected, pairwise, coverage))
        overall_rows.append(overall)
    overall = pd.DataFrame(overall_rows)
    pairwise = pd.concat([table for _, table, _ in comparisons], ignore_index=True)
    for table in [overall, pairwise]:
        available = np.isfinite(table["p-value"])
        adjusted = multipletests(table["p-value"].fillna(1), method="holm")[1]
        table["Holm-adjusted p-value"] = np.where(available, adjusted, np.nan)
    pairwise["Significant"] = pairwise["Holm-adjusted p-value"].lt(.05)
    comparisons = [(selected, pairwise.loc[pairwise["Dataset coverage"].eq(
        coverage if coverage is not None else "All")].copy(), coverage)
        for selected, _, coverage in comparisons]
    return comparisons, overall


def coverage_mean_statistics(by_dataset, combined):
    """Compare coverage means, with one averaged score per lipid per panel."""
    overall, pairwise = [], []
    for dataset in datasets + ["Combined"]:
        frame = combined if dataset == "Combined" else by_dataset.loc[by_dataset["Dataset"].eq(dataset)]
        selected = frame.loc[np.isfinite(frame[SCORE])]
        if selected["Lipid ID"].duplicated().any():
            raise ValueError("Coverage tests require one score per lipid per panel")
        counts = selected.groupby("Dataset coverage").size()
        levels = [level for level in range(1, 6) if counts.get(level, 0) >= 2]
        modeled = selected.loc[selected["Dataset coverage"].isin(levels)]
        n = len(modeled)
        result = {"Dataset": dataset, "Lipids": n, "Coverage groups tested": "; ".join(map(str, levels)),
                  "F statistic": np.nan, "Numerator df": np.nan, "Denominator df": n - 1,
                  "p-value": np.nan, "Status": "Unavailable: insufficient data or singular covariance"}
        fit = None
        if len(levels) > 1 and modeled[SCORE].nunique() > 1:
            design = np.column_stack([np.ones(n)] + [
                modeled["Dataset coverage"].eq(level).to_numpy(dtype=float) for level in levels[1:]])
            fit = OLS(modeled[SCORE].to_numpy(), design).fit().get_robustcov_results(
                cov_type="cluster", groups=modeled["Lipid ID"].to_numpy(),
                use_correction=True, df_correction=True, use_t=True)
            contrasts = np.eye(len(levels))[1:]
            covariance = contrasts @ fit.cov_params() @ contrasts.T
            if np.isfinite(covariance).all() and np.linalg.matrix_rank(covariance) == len(levels) - 1:
                test = fit.f_test(contrasts)
                result.update({"F statistic": float(np.asarray(test.fvalue).item()),
                               "Numerator df": test.df_num, "p-value": float(test.pvalue),
                               "Status": "Exploratory robust coverage mean test"})
        overall.append(result)
        for left, right in combinations(range(1, 6), 2):
            x = selected.loc[selected["Dataset coverage"].eq(left), SCORE]
            y = selected.loc[selected["Dataset coverage"].eq(right), SCORE]
            row = {"Dataset": dataset, "Coverage A": left, "Coverage B": right,
                   "Lipids A": len(x), "Lipids B": len(y), "Mean A": x.mean(), "Mean B": y.mean(),
                   "Mean difference (B - A)": np.nan, "Standard error": np.nan,
                   "t statistic": np.nan, "Denominator df": n - 1, "p-value": np.nan,
                   "Status": "Unavailable: insufficient data or singular covariance"}
            if fit is not None and left in levels and right in levels:
                contrast = np.zeros(len(levels))
                for level, sign in [(left, -1), (right, 1)]:
                    if level != levels[0]:
                        contrast[levels.index(level)] = sign
                variance = float(contrast @ fit.cov_params() @ contrast)
                if np.isfinite(variance) and variance > 0:
                    test = fit.t_test(contrast)
                    row.update({"Mean difference (B - A)": float(np.asarray(test.effect).item()),
                                "Standard error": float(np.asarray(test.sd).item()),
                                "t statistic": float(np.asarray(test.tvalue).item()),
                                "p-value": float(np.asarray(test.pvalue).item()),
                                "Status": "Exploratory robust coverage mean contrast"})
            pairwise.append(row)
    overall, pairwise = pd.DataFrame(overall), pd.DataFrame(pairwise)
    for table in [overall, pairwise]:
        available = np.isfinite(table["p-value"])
        adjusted = multipletests(table["p-value"].fillna(1), method="holm")[1]
        table["Holm-adjusted p-value"] = np.where(available, adjusted, np.nan)
    pairwise["Significant"] = pairwise["Holm-adjusted p-value"].lt(.05)
    return overall, pairwise


def median_interval(values):
    """Exploratory lipid-level percentile CI, conditional on fixed network scores."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    result = {"Lipids": len(values), "Median": np.median(values) if len(values) else np.nan,
              "Median CI lower": np.nan, "Median CI upper": np.nan,
              "CI status": "Unavailable: fewer than two scored lipids"}
    if len(values) >= 2:
        interval = bootstrap(
            (values,), np.median, confidence_level=0.95, method="percentile",
            n_resamples=BOOTSTRAPS, batch=100, random_state=np.random.default_rng(42),
        ).confidence_interval
        result.update({"Median CI lower": interval.low, "Median CI upper": interval.high,
                       "CI status": "Degenerate" if interval.low == interval.high else "Exploratory"})
    return result


def plot_median_interval(values, position, ax):
    summary = median_interval(values)
    lower, upper = summary["Median CI lower"], summary["Median CI upper"]
    if np.isfinite(lower) and np.isfinite(upper):
        # Draw endpoints directly: a percentile interval need not contain the estimate.
        ax.vlines(position, lower, upper, color="black", linewidth=1.2, zorder=5)
        ax.hlines([lower, upper], position - 0.07, position + 0.07,
                  color="black", linewidth=1.2, zorder=5)
    if np.isfinite(summary["Median"]):
        ax.plot(position, summary["Median"], "o", color="black", markersize=3, zorder=6)


def summarize_median_intervals(combined, by_dataset):
    """Summarize coverage-specific medians and pooled analyzer medians."""
    rows = []
    groups = [(dataset, by_dataset.loc[by_dataset["Dataset"].eq(dataset)])
              for dataset in datasets] + [("Combined", combined)]
    for dataset, frame in groups:
        levels = list(range(1, len(datasets) + 1))
        if dataset != "Combined":
            levels.append(None)
        for coverage in levels:
            values = (frame[SCORE] if coverage is None else
                      frame.loc[frame["Dataset coverage"].eq(coverage), SCORE]).to_numpy()
            rows.append({"Dataset": dataset, "Dataset coverage": coverage if coverage is not None else "All",
                         **median_interval(values), "Confidence level": 0.95,
                         "CI method": "Lipid-level percentile bootstrap",
                         "Bootstrap resamples": BOOTSTRAPS, "Seed": 42})
    return pd.DataFrame(rows)


def plot_comparison(paired, ax, title, bracket_ax, overall, pairwise, bracket_rows):
    scored = paired.loc[paired[SCORE].notna()]
    rng = np.random.default_rng(42)
    for coverage in range(1, len(datasets) + 1):
        group = scored.loc[scored["Dataset coverage"].eq(coverage)]
        values = group[SCORE].to_numpy()
        if len(values):
            if len(values) > 1 and np.ptp(values) > 0:
                violin = ax.violinplot([values], positions=[coverage], widths=0.7,
                                      showmedians=True, showextrema=False)
                for body in violin["bodies"]:
                    body.set_facecolor(coverage_colors[coverage])
                    body.set_edgecolor(coverage_colors[coverage])
                    body.set_alpha(0.5)
                violin["cmedians"].set_color("black")
                violin["cmedians"].set_linewidth(1)
            else:
                # Singleton and constant groups have no estimable density.
                ax.hlines(values[0], coverage - 0.2, coverage + 0.2,
                          color="black", linewidth=1)
            ax.scatter(coverage + rng.uniform(-0.12, 0.12, len(values)), values,
                       s=7, color="#347FA8",
                       alpha=0.5, linewidths=0, zorder=3)
            plot_median_interval(values, coverage, ax)
        ax.text(coverage, 1.02, f"n = {len(values)}", ha="center", va="bottom",
                transform=ax.get_xaxis_transform())
    bracket_ax.set_title(title, loc="left", pad=4, fontweight="bold")
    pvalue = overall["Holm-adjusted p-value"]
    label = f"F test: p\u2090={pvalue:.3g}" if np.isfinite(pvalue) else "F test: p\u2090=N/A"
    bracket_ax.set_title(label, loc="right", pad=4, fontsize=7)
    significant = pairwise.loc[pairwise["Significant"]].copy()
    significant["Span"] = significant["Coverage B"] - significant["Coverage A"]
    for index, (_, pair) in enumerate(significant.sort_values(["Span", "Coverage A"]).iterrows()):
        left, right = pair["Coverage A"], pair["Coverage B"]
        height = index + .25
        bracket_ax.plot([left, left, right, right], [height, height+.15, height+.15, height],
                        color="black", linewidth=.7)
        p = pair["Holm-adjusted p-value"]
        stars = "****" if p < .0001 else "***" if p < .001 else "**" if p < .01 else "*"
        bracket_ax.text((left+right)/2, height+.15, stars, ha="center", va="center", fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": .3})
    if significant.empty:
        label = "No Significant Pairs" if pairwise["p-value"].notna().any() else "Pairwise Tests Unavailable"
        bracket_ax.text(.5, .5, label, transform=bracket_ax.transAxes, ha="center", va="center", fontsize=7)
    bracket_ax.set_xlim(.5, 5.5)
    bracket_ax.set_ylim(0, bracket_rows+.3)
    bracket_ax.set_axis_off()
    ax.set_xticks(range(1, len(datasets) + 1))
    ax.set_xlim(0.5, len(datasets) + 0.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1], labels=["0", "0.25", "0.5", "0.75", "1"])
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def plot_dataset_comparison(paired, combined, overall, pairwise):
    fig = plt.figure(figsize=(6.5, 6.75), layout="constrained")
    grid = fig.add_gridspec(3, 2)
    significant = pairwise.loc[pairwise["Significant"]]
    bracket_rows = max(1, int(significant.groupby("Dataset").size().max())) if len(significant) else 1
    axes = []
    for index, dataset in enumerate(["Combined"] + datasets):
        cell = grid[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=[.12 + .14 * bracket_rows, 1.55], hspace=.02)
        bracket_ax = fig.add_subplot(cell[0])
        ax = fig.add_subplot(cell[1], sharex=axes[0] if axes else None, sharey=axes[0] if axes else None)
        axes.append(ax)
        group = combined if dataset == "Combined" else paired.loc[paired["Dataset"].eq(dataset)]
        plot_comparison(group, ax, f"{chr(65 + index)}  {DATASET_LABELS.get(dataset, dataset)}",
                        bracket_ax, overall.loc[overall["Dataset"].eq(dataset)].iloc[0],
                        pairwise.loc[pairwise["Dataset"].eq(dataset)], bracket_rows)
        ax.tick_params(labelbottom=index >= 4, labelleft=index % 2 == 0)
    fig.supxlabel("Number Of Datasets Confirming Lipid", fontsize=9)
    fig.supylabel("Mean Normalized Similarity Score", fontsize=9)
    return fig


def plot_coverage_comparison(selected, statistics, coverage, ax, overall):
    rng = np.random.default_rng(42)
    color = coverage_colors[coverage] if coverage is not None else "#347FA8"
    for position, dataset in enumerate(datasets):
        values = selected.loc[selected["Dataset"].eq(dataset), SCORE].to_numpy()
        if len(values):
            if len(values) > 1 and np.ptp(values) > 0:
                violin = ax.violinplot([values], positions=[position], widths=0.7,
                                      showmedians=True, showextrema=False)
                for body in violin["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(color)
                    body.set_alpha(0.5)
                violin["cmedians"].set_color("black")
                violin["cmedians"].set_linewidth(1)
            else:
                ax.hlines(values[0], position - 0.2, position + 0.2,
                          color="black", linewidth=1)
            ax.scatter(position + rng.uniform(-0.12, 0.12, len(values)), values,
                       s=7, color="#347FA8", alpha=0.5, linewidths=0, zorder=3)
            plot_median_interval(values, position, ax)
        ax.text(position, 1.04, f"n = {len(values)}", ha="center", va="bottom")
    significant = statistics.loc[statistics["Significant"]].copy()
    significant["Span"] = [datasets.index(row["Dataset"])
                           - datasets.index(row["Reference Dataset"])
                           for _, row in significant.iterrows()]
    for level, (_, row) in enumerate(significant.sort_values("Span").iterrows()):
        left = datasets.index(row["Reference Dataset"])
        right = datasets.index(row["Dataset"])
        height = 1.17 + 0.10 * level
        ax.plot([left, left, right, right], [height - 0.025, height, height, height - 0.025],
                color="black", linewidth=0.7)
        p = row["Holm-adjusted p-value"]
        stars = "****" if p < 0.0001 else "***" if p < 0.001 else "**" if p < 0.01 else "*"
        ax.text((left + right) / 2, height, stars, ha="center", va="center", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.3})
    title = ("Combined" if coverage is None else
             f"{coverage} Dataset" + ("s" if coverage > 1 else ""))
    ax.set_title(title, loc="left", pad=8, fontweight="bold")
    pvalue = overall["Holm-adjusted p-value"]
    label = f"F test: p\u2090={pvalue:.3g}" if np.isfinite(pvalue) else "F test: p\u2090=N/A"
    ax.set_title(label, loc="right", pad=8, fontsize=7)
    ax.set_xticks(range(len(datasets)), labels=[DATASET_LABELS[d] for d in datasets],
                  rotation=35, ha="right")
    ax.set_xlim(-0.5, len(datasets) - 0.5)
    ax.set_ylim(-0.03, 1.13 + 0.10 * len(significant))
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1], labels=["0", "0.25", "0.5", "0.75", "1"])
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def plot_coverage_panels(comparisons, overall):
    fig, axes = plt.subplots(3, 2, figsize=(6.5, 6.75), sharex=True, sharey=True,
                             layout="constrained")
    upper = 1.13 + 0.10 * max(int(statistics["Significant"].sum())
                             for _, statistics, _ in comparisons)
    for index, ((selected, statistics, coverage), ax) in enumerate(zip(comparisons, axes.flat)):
        plot_coverage_comparison(selected, statistics, coverage, ax, overall.iloc[index])
        ax.set_title(f"{chr(65 + index)}  {ax.get_title(loc='left')}",
                     loc="left", pad=8, fontweight="bold")
        ax.set_ylim(-0.03, upper)
        # Show analyzer labels only on the bottom row.
        ax.tick_params(axis="x", labelbottom=index >= len(comparisons) - 2)
    fig.supxlabel("Analyzer", fontsize=9)
    fig.supylabel("Mean Normalized Similarity Score", fontsize=9)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    source = (pd.read_csv(args.input) if args.input.suffix.casefold() == ".csv"
              else pd.read_excel(args.input, sheet_name="Sheet1"))
    scored, _ = recompute_adduct_networks(source)
    points = prepare_strength(scored, all_ids=False)
    paired = prepare_comparison(source, points)
    by_dataset = prepare_dataset_comparison(paired, points)
    coverage_overall, coverage_pairwise = coverage_mean_statistics(by_dataset, paired)
    dataset_fig = plot_dataset_comparison(by_dataset, paired, coverage_overall, coverage_pairwise)
    output_dir = args.output_dir or args.input.parent / "Final Figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        dataset_fig.savefig(output_dir / f"SI_6.{extension}", dpi=600)
    comparisons, overall = analyzer_comparisons(by_dataset)
    for selected, statistics, coverage in comparisons:
        print(statistics.to_string(index=False))
    comparison_fig = plot_coverage_panels(comparisons, overall)
    for extension in ["png", "pdf"]:
        comparison_fig.savefig(output_dir / f"SI_5.{extension}", dpi=600)
    print(overall.to_string(index=False))
    print("Exploratory robust ANOVA-style F tests of means and two-sided t-tests of mean contrasts, clustered by lipid. "
          "Holm correction across 6 overall tests and across 60 pairwise tests. "
          "Stars: * <0.05, ** <0.01, *** <0.001, **** <0.0001 (adjusted p).")
    print("Black dots/bars: medians and pointwise 95% lipid-level percentile-bootstrap CIs "
          "(9,999 resamples, seed 42), conditional on fixed networks; not run-level CIs.")
    print(paired.groupby("Dataset coverage")[SCORE].agg(["count", "mean", "median"]).to_string())
    print(f"Confirmed lipids without scores: {paired[SCORE].isna().sum()}")
    if args.no_show:
        plt.close(dataset_fig)
        plt.close(comparison_fig)
    else:
        plt.show(block=True)


if __name__ == "__main__":
    main()
