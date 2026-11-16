"""Violin plots of normalized MS1/MS2 intensity versus confirmed dataset coverage.

Accepted, non-provisional IDs only; original ontology classes remain separate.
MS1 intensity is the nearest precursor peak in MS1 isotopic spectrum, using
0.05 m/z for OTOT/OTIT and 0.3 for ITIT/Stellar. MS2 intensity is the sum of
exported MS/MS intensities, with no threshold or precursor removal. No score,
network construction, or network-size restriction is used.

Collapse duplicate features to median intensity per dataset/original-class/
polarity/adduct/lipid identity (name before '|', ignoring case and whitespace).
Normalize each median to 100 * intensity / maximum within that SAME dataset/
class/polarity/adduct, separately for MS1 and MS2, BEFORE coverage aggregation.
All-zero groups retain zero; missing/invalid intensities are omitted, never
zero-filled. Coverage counts distinct datasets confirming each original-class/
lipid identity across adducts/polarities, independently of intensity availability.

Average normalized adduct/polarity values for each lipid WITHIN each dataset.
Never average across datasets. For each of MS1/MS2, show five dataset panels
plus a Combined panel pooling all dataset-specific lipid values, displayed first.
A lipid seen
in several datasets contributes once per dataset to Combined. Counts in that
panel are lipid-dataset observations, not unique lipid identities. X is the
number of confirming datasets (1-5); Y is within-dataset mean normalized
intensity (%). Coverage is still measured across all five datasets. This is a
relative within-class/adduct/polarity comparison, not absolute abundance.
Violins show distributions/medians; jittered dots show every lipid. Constant
and singleton groups use a horizontal mark instead of a density estimate.

Run directly; set MPLBACKEND=Agg for headless export. Save PNG (600 dpi), editable
PDF figures only, named SI_7_MS1 and SI_7_MS2.
The main figures use a shared linear percentage axis. Optional log-y exports
are controlled by EXPORT_LOG_Y. Log violins estimate density in displayed log
coordinates; medians/dots retain their original values. If zeros occur, use a
shared symmetric-log scale with a linear region below the smallest positive
intensity, retaining zeros without adding a pseudocount.

Statistics use original percentage intensities, identically for both displays.
OLS tests mean intensity differences across coverage categories, with CR1
covariance clustered by lipid identity and G-1 denominator degrees of freedom.
The overall ANOVA-style F test jointly compares means, not variances;
the pairwise two-sided t contrasts compare individual mean differences.
Combined includes dataset fixed effects and keeps all lipid-dataset observations.
Report the overall coverage F test (Holm-adjusted across 12 planned MS1/MS2 panels)
on each panel and calculate all 120 planned two-sided pairwise t contrasts with
one Holm correction across both intensity types and all panels. Missing or
singleton coverage groups have unavailable contrasts, counted as p=1 in correction.
Tests concern means; violins mark medians. Results remain exploratory: clustering
does not account for dependencies between different lipids sharing normalization
denominators, or establish independent instrument-run replication. Covariance
conditions on the observed normalization and coverage assignments.
Show brackets only for pairwise Holm-adjusted p < 0.05: * p < 0.05,
** p < 0.01, *** p < 0.001, **** p < 0.0001. Omitted pairs are either
not significant or unavailable; all planned comparisons remain in memory.
Bracket strips sit above the intensity axes and use the same x positions.
Both displays use identical tests and annotations. The compact 6.5 x 6.75-inch
canvas leaves page space for a caption while retaining separate bracket strips.
"""

import os
from pathlib import Path

import matplotlib
if __name__ == "__main__":
    matplotlib.use(os.environ.get("MPLBACKEND", "WxAgg"))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations
from statsmodels.regression.linear_model import OLS

INPUT_FILE = Path(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline\Combined_Fragment_Assignment_Colored.xlsx")
OUTPUT_DIR = INPUT_FILE.parent / "Final Figures"
BASENAME = "SI_7"
DATASETS = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
PANEL_ORDER = ["Combined"] + DATASETS
EXPORT_LOG_Y = False
DATASET_LABELS = {"OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT",
                  "Stellar": "Stellar", "Stellarboth": "Stellar (\u00b1)"}
MS1_TOLERANCES = {"OTOT": 0.05, "OTIT": 0.05, "ITIT": 0.3,
                  "Stellar": 0.3, "Stellarboth": 0.3}
GROUPS = ["Dataset", "Ontology", "Polarity", "Adduct type"]
RELATIVE = "Normalized intensity (%)"
RAW = "Intensity"
COLORS = ["#D9EAF4", "#B8D8E8", "#8FC1DA", "#5FA6C8", "#347FA8"]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def spectrum_intensity(value, precursor=None, tolerance=None):
    """Return (intensity, status); nearest precursor peak or full spectrum sum."""
    if not isinstance(value, str) or not value.strip():
        return np.nan, "Missing spectrum"
    try:
        peaks = np.asarray([token.split(":") for token in value.split()], dtype=float)
    except (ValueError, TypeError):
        return np.nan, "Malformed spectrum"
    if (peaks.ndim != 2 or peaks.shape[1] != 2 or not np.isfinite(peaks).all()
            or (peaks[:, 0] <= 0).any() or (peaks[:, 1] < 0).any()):
        return np.nan, "Invalid m/z or intensity"
    if tolerance is not None:
        if not np.isfinite(tolerance) or tolerance <= 0:
            return np.nan, "Invalid MS1 tolerance"
        if not np.isfinite(precursor) or precursor <= 0:
            return np.nan, "Invalid precursor m/z"
        distances = np.abs(peaks[:, 0] - precursor)
        nearest = int(np.argmin(distances))
        if distances[nearest] > tolerance:
            return np.nan, "No MS1 peak within tolerance"
        return float(peaks[nearest, 1]), "OK"
    with np.errstate(over="ignore"):
        total = float(peaks[:, 1].sum())
    return (total, "OK") if np.isfinite(total) else (np.nan, "Intensity sum overflow")


def prepare_tables(source):
    """Normalize intensities before joining intensity-independent ID coverage."""
    required = set(GROUPS + ["Metabolite name", "Overall Match Status",
                             "Average Mz", "MS1 isotopic spectrum", "MS/MS spectrum"])
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))
    frame = source.copy().reset_index(drop=True)
    frame["Input row (zero-based)"] = frame.index
    for column in GROUPS + ["Metabolite name"]:
        frame[column] = frame[column].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.split().str.join(" ")
    frame["Dataset"] = frame["Dataset"].str.casefold().map({d.casefold(): d for d in DATASETS})
    frame["Polarity"] = frame["Polarity"].str.casefold()
    canonical = {name.casefold(): name for name in sorted(frame["Ontology"].unique(), reverse=True)}
    frame["Ontology"] = frame["Ontology"].str.casefold().map(canonical)
    frame["Lipid Name"] = frame["Metabolite name"].str.split("|", n=1).str[0].str.strip()
    frame["Lipid ID"] = frame["Ontology"].str.casefold() + "|" + frame["Lipid Name"].str.casefold()
    frame["Average Mz"] = pd.to_numeric(frame["Average Mz"], errors="coerce")
    accepted = frame["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    if "Validation Status" in frame:
        accepted &= frame["Validation Status"].fillna("").astype(str).str.replace(
            "\u00a0", " ", regex=False).str.strip().str.casefold().ne("provisional")
    identified = (accepted & frame["Dataset"].notna() & frame["Ontology"].ne("")
                  & frame["Lipid Name"].ne(""))
    coverage = frame.loc[identified].groupby("Lipid ID", sort=True).agg(
        **{"Lipid Name": ("Lipid Name", "first"), "Ontology": ("Ontology", "first"),
           "Dataset coverage": ("Dataset", "nunique")}).reset_index()
    eligible = identified & frame["Polarity"].isin(["positive", "negative"]) & frame["Adduct type"].ne("")
    frame["Eligible identification"] = eligible
    frame["MS1 tolerance (m/z)"] = frame["Dataset"].map(MS1_TOLERANCES)
    ms1 = [spectrum_intensity(value, mz, tol) for value, mz, tol in zip(
        frame["MS1 isotopic spectrum"], frame["Average Mz"], frame["MS1 tolerance (m/z)"])]
    ms2 = [spectrum_intensity(value) for value in frame["MS/MS spectrum"]]
    tables = []
    for kind, calculated in [("MS1", ms1), ("MS2", ms2)]:
        frame[f"{kind} intensity"] = [item[0] for item in calculated]
        frame[f"{kind} status"] = [item[1] for item in calculated]
        frame[f"{kind} included"] = eligible & frame[f"{kind} status"].eq("OK")
        selected = frame.loc[frame[f"{kind} included"]]
        points = selected.groupby(GROUPS + ["Lipid ID"], sort=True).agg(
            **{"Lipid Name": ("Lipid Name", "first"), RAW: (f"{kind} intensity", "median"),
               "Feature rows": (f"{kind} intensity", "size")}).reset_index()
        points["Intensity type"] = kind
        points["Normalization maximum"] = points.groupby(GROUPS)[RAW].transform("max")
        maximum = points["Normalization maximum"]
        points[RELATIVE] = points[RAW] / maximum.where(maximum.gt(0)) * 100
        points.loc[maximum.eq(0), RELATIVE] = 0.0
        tables.append(points)
    points = pd.concat(tables, ignore_index=True)
    if points.empty:
        raise ValueError("No eligible intensities")
    audit = points.groupby(["Intensity type"] + GROUPS, sort=True).agg(
        **{"Lipids": ("Lipid ID", "size"), "Maximum intensity": (RAW, "max"),
           "Normalized maximum (%)": (RELATIVE, "max")}).reset_index()
    audit["Expected maximum (%)"] = np.where(audit["Maximum intensity"].gt(0), 100., 0.)
    audit["Pass"] = np.isclose(audit["Normalized maximum (%)"], audit["Expected maximum (%)"])
    if not audit["Pass"].all():
        raise ValueError("Intensity normalization audit failed")
    by_dataset = points.groupby(["Intensity type", "Lipid ID", "Dataset"], sort=True)[RELATIVE].mean().reset_index()
    by_dataset = by_dataset.merge(coverage, on="Lipid ID", validate="many_to_one")
    # Keep every dataset-specific lipid observation in the pooled panel/export.
    combined = by_dataset.copy()
    return points, combined, by_dataset, coverage, frame, audit


def summarize_coverage(combined):
    rows = []
    for kind in ["MS1", "MS2"]:
        population = combined.loc[combined["Intensity type"].eq(kind)]
        for dataset in PANEL_ORDER:
            selected = population if dataset == "Combined" else population.loc[population["Dataset"].eq(dataset)]
            for count in range(1, 6):
                group = selected.loc[selected["Dataset coverage"].eq(count)]
                values = group[RELATIVE]
                rows.append({"Intensity type": kind, "Dataset": dataset,
                             "Dataset coverage": count, "Lipid-dataset observations": len(values),
                             "Unique lipids": group["Lipid ID"].nunique(),
                             "Mean normalized intensity (%)": values.mean(),
                             "Median normalized intensity (%)": values.median(),
                             "Q1 (%)": values.quantile(.25), "Q3 (%)": values.quantile(.75)})
    return pd.DataFrame(rows)


def holm_adjust(table, column):
    """Include unavailable planned tests as p=1 in the correction family."""
    values = table[column].to_numpy(dtype=float)
    planned = np.where(np.isfinite(values), values, 1.)
    order = np.argsort(planned)
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(1., np.maximum.accumulate(
        planned[order] * np.arange(len(values), 0, -1)))
    adjusted[~np.isfinite(values)] = np.nan
    return adjusted


def coverage_statistics(combined):
    """Coverage mean tests with lipid-clustered covariance; dataset-adjust pooled."""
    overall, pairwise = [], []
    for kind in ["MS1", "MS2"]:
        population = combined.loc[combined["Intensity type"].eq(kind)]
        for dataset in PANEL_ORDER:
            selected = population if dataset == "Combined" else population.loc[population["Dataset"].eq(dataset)]
            counts = selected.groupby("Dataset coverage").size()
            levels = sorted(counts.index)
            clusters = selected["Lipid ID"].nunique()
            row = {"Intensity type": kind, "Dataset": dataset, "Observations": len(selected),
                   "Lipid clusters": clusters,
                   "F statistic": np.nan, "Numerator df": np.nan, "Denominator df": clusters - 1,
                   "p-value": np.nan, "Status": "Unavailable: insufficient data or singular model"}
            fit = None
            columns = []
            # Singleton strata have no useful variance estimate for their contrasts.
            if (len(levels) > 1 and clusters > 1 and selected[RELATIVE].nunique() > 1
                    and all(counts[level] >= 2 for level in levels)):
                design = np.ones((len(selected), 1))
                columns = levels[1:]
                for level in columns:
                    design = np.column_stack((design, selected["Dataset coverage"].eq(level).astype(float)))
                if dataset == "Combined":
                    for name in sorted(selected["Dataset"].unique())[1:]:
                        design = np.column_stack((design, selected["Dataset"].eq(name).astype(float)))
                if len(selected) > design.shape[1] and np.linalg.matrix_rank(design) == design.shape[1]:
                    model = OLS(selected[RELATIVE].to_numpy(), design).fit()
                    fit = model.get_robustcov_results(
                        cov_type="cluster", groups=selected["Lipid ID"].to_numpy(),
                        use_correction=True, df_correction=True, use_t=True)
                    contrasts = np.eye(design.shape[1])[1:len(columns) + 1]
                    covariance = contrasts @ fit.cov_params() @ contrasts.T
                    if np.isfinite(covariance).all() and np.linalg.matrix_rank(covariance) == len(columns):
                        test = fit.f_test(contrasts)
                        row.update({"F statistic": float(np.asarray(test.fvalue).item()),
                                    "Numerator df": test.df_num, "p-value": float(test.pvalue),
                                    "Status": "Exploratory clustered mean test"})
            overall.append(row)
            for left, right in combinations(range(1, 6), 2):
                comparison = {"Intensity type": kind, "Dataset": dataset,
                              "Coverage A": left, "Coverage B": right,
                              "Observations A": int(counts.get(left, 0)),
                              "Observations B": int(counts.get(right, 0)),
                              "Mean difference (B - A)": np.nan, "p-value": np.nan,
                              "Standard error": np.nan, "t statistic": np.nan,
                              "Denominator df": clusters - 1,
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
                                           "Standard error": float(np.asarray(test.sd).item()),
                                           "t statistic": float(np.asarray(test.tvalue).item()),
                                           "p-value": float(np.asarray(test.pvalue).item()),
                                           "Status": "Exploratory clustered mean contrast"})
                pairwise.append(comparison)
    overall, pairwise = pd.DataFrame(overall), pd.DataFrame(pairwise)
    for table in [overall, pairwise]:
        table["Holm-adjusted p-value"] = holm_adjust(table, "p-value")
    return overall, pairwise


def significance_label(pvalue):
    """Return stars for a finite, significant Holm-adjusted p-value."""
    if not np.isfinite(pvalue) or pvalue >= .05:
        return ""
    return "****" if pvalue < .0001 else "***" if pvalue < .001 else "**" if pvalue < .01 else "*"


def plot_comparison(combined, kind, statistics, pairwise, log_y=False):
    comparisons = pairwise.loc[pairwise["Intensity type"].eq(kind)].copy()
    comparisons["Significance"] = comparisons["Holm-adjusted p-value"].map(significance_label)
    significant = comparisons.loc[comparisons["Significance"].ne("")]
    bracket_rows = max(1, int(significant.groupby("Dataset").size().max())) if len(significant) else 1
    fig = plt.figure(figsize=(6.5, 6.75), layout="constrained")
    fig.get_layout_engine().set(h_pad=.025, w_pad=.025, hspace=.025, wspace=.025)
    grid = fig.add_gridspec(3, 2)
    axes, bracket_axes = [], []
    for index in range(6):
        cell = grid[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=[.12 + .14 * bracket_rows, 1.55], hspace=.02)
        bracket_axes.append(fig.add_subplot(cell[0]))
        ax = fig.add_subplot(cell[1], sharex=axes[0] if axes else None,
                             sharey=axes[0] if axes else None)
        ax.tick_params(labelbottom=index >= 4, labelleft=index % 2 == 0)
        axes.append(ax)
    rng = np.random.default_rng(42)
    population = combined.loc[combined["Intensity type"].eq(kind)]
    positive = population.loc[population[RELATIVE].gt(0), RELATIVE]
    minimum = float(positive.min()) if len(positive) else 1.
    for index, (dataset, ax, bracket_ax) in enumerate(zip(PANEL_ORDER, axes, bracket_axes)):
        if log_y:
            if population[RELATIVE].eq(0).any():
                ax.set_yscale("symlog", linthresh=minimum)
                ax.set_ylim(-minimum * .1, 200)
            else:
                ax.set_yscale("log")
                ax.set_ylim(minimum / 2, 200)
        selected = population if dataset == "Combined" else population.loc[population["Dataset"].eq(dataset)]
        for count, color in zip(range(1, 6), COLORS):
            values = selected.loc[selected["Dataset coverage"].eq(count), RELATIVE].to_numpy()
            if len(values):
                if len(values) > 1 and np.ptp(values) > 0:
                    transform = ax.yaxis.get_transform()
                    displayed = transform.transform(values) if log_y else values
                    violin = ax.violinplot([displayed], positions=[count], widths=.7,
                                          showmedians=True, showextrema=False)
                    for body in violin["bodies"]:
                        if log_y:
                            for path in body.get_paths():
                                path.vertices[:, 1] = transform.inverted().transform(path.vertices[:, 1])
                        body.set_facecolor(color)
                        body.set_edgecolor(color)
                        body.set_alpha(.6)
                    violin["cmedians"].set_color("black")
                    if log_y:
                        segments = violin["cmedians"].get_segments()
                        for segment in segments:
                            segment[:, 1] = np.median(values)
                        violin["cmedians"].set_segments(segments)
                else:
                    ax.hlines(values[0], count - .2, count + .2, color="black", linewidth=1)
                ax.scatter(count + rng.uniform(-.12, .12, len(values)), values,
                           s=5, color="#347FA8", alpha=.35, linewidths=0, zorder=3)
            ax.text(count, .96, f"n={len(values)}", ha="center", va="center", fontsize=7,
                    transform=ax.get_xaxis_transform())
        label = DATASET_LABELS.get(dataset, dataset)
        bracket_ax.set_title(f"{chr(65 + index)}  {label}", loc="left", fontweight="bold")
        result = statistics.loc[statistics["Intensity type"].eq(kind) & statistics["Dataset"].eq(dataset)].iloc[0]
        pvalue = result["Holm-adjusted p-value"]
        p_text = f"F test: p\u2090={pvalue:.2g}" if np.isfinite(pvalue) else "F test: p\u2090=N/A"
        bracket_ax.set_title(p_text, loc="right", fontsize=7)
        pairs = significant.loc[significant["Dataset"].eq(dataset)].copy()
        pairs["Span"] = pairs["Coverage B"] - pairs["Coverage A"]
        pairs = pairs.sort_values(["Span", "Coverage A", "Coverage B"])
        for row, (_, pair) in enumerate(pairs.iterrows()):
            left, right = pair["Coverage A"], pair["Coverage B"]
            height = row + .25
            bracket_ax.plot([left, left, right, right],
                            [height, height + .15, height + .15, height],
                            color="black", linewidth=.7)
            bracket_ax.text((left + right) / 2, height + .15, pair["Significance"],
                            ha="center", va="center", fontsize=8,
                            bbox={"facecolor": "white", "edgecolor": "none", "pad": .3})
        if pairs.empty:
            available = comparisons.loc[comparisons["Dataset"].eq(dataset), "Holm-adjusted p-value"].notna().any()
            bracket_ax.text(.5, .5, "No Significant Pairs" if available else "Pairwise Tests Unavailable",
                            transform=bracket_ax.transAxes, ha="center", va="center", fontsize=7)
        bracket_ax.set_xlim(.5, 5.5)
        bracket_ax.set_ylim(0, bracket_rows + .3)
        bracket_ax.set_axis_off()
        ax.set_xlim(.5, 5.5)
        if not log_y:
            ax.set_ylim(0, 112)
            ax.set_yticks([0, 50, 100])
        ax.set_xticks(range(1, 6))
        ax.spines[["top", "right"]].set_visible(False)
    fig.supxlabel("Number Of Datasets Confirming Lipid", fontsize=9)
    quantity = "Summed MS2" if kind == "MS2" else "MS1 Precursor"
    scale = "; Log Scale" if log_y else ""
    fig.supylabel(f"Mean Normalized {quantity} Intensity (%{scale})", fontsize=9)
    return fig


if __name__ == "__main__":
    source = pd.read_excel(INPUT_FILE, sheet_name="Sheet1")
    points, combined, by_dataset, coverage, features, normalization = prepare_tables(source)
    statistics, pairwise = coverage_statistics(combined)
    pairwise["Significance"] = pairwise["Holm-adjusted p-value"].map(significance_label)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figures = []
    for kind in ["MS1", "MS2"]:
        for log_y in ([False, True] if EXPORT_LOG_Y else [False]):
            figure = plot_comparison(combined, kind, statistics, pairwise, log_y=log_y)
            figures.append(figure)
            suffix = "_LogY" if log_y else ""
            for extension in ["png", "pdf"]:
                figure.savefig(OUTPUT_DIR / f"{BASENAME}_{kind}{suffix}.{extension}", dpi=600)
    print(summarize_coverage(combined).to_string(index=False))
    print(f"Saved outputs to {OUTPUT_DIR}")
    if matplotlib.get_backend().casefold() == "agg":
        for figure in figures:
            plt.close(figure)
    else:
        plt.show(block=True)
