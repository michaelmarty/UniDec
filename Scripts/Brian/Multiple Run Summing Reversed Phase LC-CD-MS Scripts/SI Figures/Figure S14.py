"""Create the Rabbit IgG delta plots used in Figure S14.

Each plotted value is the replicate mean at an injection count minus the
replicate mean at injection 21. Error bars are SEM in the original metric
units. Peak-area fractions are normalized over the retained target peaks in
each UniDec peak-parameter file.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


REPLICATE_FOLDERS = [
    Path(
        r"Z:\Group Share\BHT\Temp\Paper\Stats\High Throughput"
        r"\Amgen Rabbit IgG Rep 1"
    ),
    Path(
        r"Z:\Group Share\BHT\Temp\Paper\Stats\High Throughput"
        r"\Amgen Rabbit IgG Rep 2"
    ),
    Path(
        r"Z:\Group Share\BHT\Temp\Paper\Stats\High Throughput"
        r"\Amgen Rabbit IgG Rep 3"
    ),
]

TARGET_MASSES = {
    "142,631 Da": 142631.0,
    "142,828 Da": 142828.0,
    "143,028 Da": 143028.0,
    "143,286 Da": 143286.0,
    "143,584 Da": 143584.0,
    "143,752 Da": 143752.0,
    "143,920 Da": 143920.0,
    "144,027 Da": 144027.0,
    "144,243 Da": 144243.0,
    "144,536 Da": 144536.0,
    "144,704 Da": 144704.0,
    "144,864 Da": 144864.0,
    "145,017 Da": 145017.0,
}

MASS_TOLERANCE_DA = 100.0
EXPECTED_INJECTIONS = tuple(range(1, 22))
REFERENCE_INJECTION = 21
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S14.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S14.pdf"
OUTPUT_CSV = OUTPUT_DIRECTORY / "Figure S14_results.csv"

BASE_FONT_SIZE = 8.5
PANEL_TITLE_FONT_SIZE = 10
SUBPLOT_LABEL_FONT_SIZE = 11
LEGEND_FONT_SIZE = 6
FIGURE_SIZE_INCHES = (7.0, 4.0)
AXIS_STROKE_WIDTH = 0.8

REQUIRED_COLUMNS = {
    "Mass",
    "Area",
    "MassCentroid",
    "MassFWHM",
}

METRICS = {
    "centroid": {
        "title": "Mass Centroid",
        "source": "Mass Centroid (Da)",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION} (Da)",
        "minimum_limit": 1.0,
    },
    "fwhm": {
        "title": "FWHM",
        "source": "FWHM (Da)",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION} (Da)",
        "minimum_limit": 1.0,
    },
    "area_fraction": {
        "title": "Sum-Normalized\nPeak Area Fraction",
        "source": "Sum-Normalized Peak Area Fraction",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION}",
        "minimum_limit": 0.01,
    },
}


def extract_injection_number(path: Path) -> int:
    """Extract an injection count from a ``merged_N_peakparam.dat`` name."""
    match = re.search(r"merged_(\d+)_", path.name)
    if match is None:
        raise ValueError(f"Could not determine injection number from: {path.name}")
    return int(match.group(1))


def load_peak_parameters(path: Path) -> pd.DataFrame:
    """Load and validate one UniDec peak-parameter file."""
    table = pd.read_csv(path, sep=r"\s+", encoding="utf-8")
    missing_columns = REQUIRED_COLUMNS.difference(table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing columns ({missing}) in: {path}")

    table = table.copy()
    for column in REQUIRED_COLUMNS:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if table[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError(f"Non-numeric required values found in: {path}")
    return table


def match_target_peaks(
    table: pd.DataFrame,
    path: Path,
) -> dict[str, pd.Series]:
    """Return the closest in-tolerance row for each available target peak."""
    matched: dict[str, pd.Series] = {}
    for peak_name, target_mass in TARGET_MASSES.items():
        mass_error = (table["Mass"] - target_mass).abs()
        candidates = table.loc[mass_error <= MASS_TOLERANCE_DA]
        if candidates.empty:
            print(
                f"Warning: No {peak_name} peak within {MASS_TOLERANCE_DA:g} Da "
                f"in {path.name}"
            )
            continue
        closest_index = mass_error.loc[candidates.index].idxmin()
        matched[peak_name] = table.loc[closest_index]
    return matched


def collect_measurements() -> pd.DataFrame:
    """Collect matched Rabbit IgG peak measurements from every replicate."""
    records: list[dict[str, float | int | str]] = []

    for replicate, folder in enumerate(REPLICATE_FOLDERS, start=1):
        if not folder.is_dir():
            raise FileNotFoundError(f"Replicate folder not found: {folder}")

        files = list(folder.glob("*_peakparam.dat"))
        files.sort(key=extract_injection_number)
        found_injections = tuple(extract_injection_number(path) for path in files)
        if found_injections != EXPECTED_INJECTIONS:
            raise ValueError(
                f"Replicate {replicate} must contain injections 1--21 exactly once; "
                f"found {found_injections}"
            )

        for path in files:
            injection_number = extract_injection_number(path)
            table = load_peak_parameters(path)
            matched_peaks = match_target_peaks(table, path)
            if not matched_peaks:
                raise ValueError(f"No target peaks were matched in: {path}")

            matched_area_sum = sum(
                float(peak["Area"]) for peak in matched_peaks.values()
            )
            if matched_area_sum <= 0.0:
                raise ValueError(f"Matched peak areas do not have a positive sum: {path}")

            for peak_name, peak in matched_peaks.items():
                peak_area = float(peak["Area"])
                if peak_area < 0.0:
                    raise ValueError(f"Negative matched peak area found in: {path}")
                records.append(
                    {
                        "Peak": peak_name,
                        "Target Mass (Da)": TARGET_MASSES[peak_name],
                        "Injection Number": injection_number,
                        "Replicate": replicate,
                        "Matched Mass (Da)": float(peak["Mass"]),
                        "Mass Centroid (Da)": float(peak["MassCentroid"]),
                        "FWHM (Da)": float(peak["MassFWHM"]),
                        "Sum-Normalized Peak Area Fraction": (
                            peak_area / matched_area_sum
                        ),
                        "Source File": str(path),
                    }
                )

    measurements = pd.DataFrame.from_records(records)
    if measurements.empty:
        raise ValueError("No Rabbit IgG target measurements were matched")
    return measurements


def standard_error(values: np.ndarray) -> float:
    """Return SEM; a single available replicate has an SEM of zero."""
    if values.size < 2:
        return 0.0
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def create_delta_results(measurements: pd.DataFrame) -> pd.DataFrame:
    """Calculate replicate means, SEM, and deltas from injection 21."""
    rows: list[dict[str, float | int | str]] = []

    for peak_name, target_mass in TARGET_MASSES.items():
        peak_data = measurements.loc[measurements["Peak"] == peak_name]
        for metric_name, metric in METRICS.items():
            reference_values = peak_data.loc[
                peak_data["Injection Number"] == REFERENCE_INJECTION,
                metric["source"],
            ].to_numpy(dtype=float)
            if reference_values.size == 0:
                print(
                    f"Warning: No injection {REFERENCE_INJECTION} reference for "
                    f"{peak_name} ({metric['title']}); series omitted"
                )
                continue
            reference_mean = float(np.mean(reference_values))

            for injection_number in EXPECTED_INJECTIONS:
                values = peak_data.loc[
                    peak_data["Injection Number"] == injection_number,
                    metric["source"],
                ].to_numpy(dtype=float)
                if values.size == 0:
                    continue
                mean = float(np.mean(values))
                rows.append(
                    {
                        "Peak": peak_name,
                        "Target Mass (Da)": target_mass,
                        "Metric": metric["title"].replace("\n", " "),
                        "Reference Injection": REFERENCE_INJECTION,
                        "Reference Mean": reference_mean,
                        "Injection Number": injection_number,
                        "Replicate Count": int(values.size),
                        "Mean": mean,
                        "Delta": mean - reference_mean,
                        "SEM": standard_error(values),
                    }
                )

    results = pd.DataFrame.from_records(rows)
    if results.empty:
        raise ValueError("No delta results could be calculated")
    return results


def get_peak_colors() -> dict[str, tuple[float, ...]]:
    """Assign each target peak a stable color from the tab20 palette."""
    colors = plt.get_cmap("tab20").colors
    return {
        peak_name: colors[index]
        for index, peak_name in enumerate(TARGET_MASSES)
    }


def set_delta_axis_limits(
    axis: plt.Axes,
    metric_results: pd.DataFrame,
    minimum_limit: float,
) -> None:
    """Set symmetric limits that contain every mean +/- SEM and zero."""
    lower = metric_results["Delta"].to_numpy(dtype=float) - metric_results[
        "SEM"
    ].to_numpy(dtype=float)
    upper = metric_results["Delta"].to_numpy(dtype=float) + metric_results[
        "SEM"
    ].to_numpy(dtype=float)
    maximum = max(
        float(np.max(np.abs(lower))),
        float(np.max(np.abs(upper))),
        minimum_limit,
    )
    axis.set_ylim(-maximum * 1.10, maximum * 1.10)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=7, symmetric=True))


def make_figure(results: pd.DataFrame) -> plt.Figure:
    """Create the three-panel Rabbit IgG delta plot."""
    figure, axes = plt.subplots(
        1,
        len(METRICS),
        squeeze=False,
        sharex=True,
        figsize=FIGURE_SIZE_INCHES,
    )
    axes = axes[0]
    peak_colors = get_peak_colors()
    legend_handles = []
    legend_labels = []

    for panel_index, (axis, (metric_name, metric)) in enumerate(
        zip(axes, METRICS.items())
    ):
        metric_title = metric["title"].replace("\n", " ")
        metric_results = results.loc[results["Metric"] == metric_title]

        axis.set_title(
            metric["title"],
            fontsize=PANEL_TITLE_FONT_SIZE,
            fontweight="bold",
            pad=14,
        )
        axis.text(
            0.97,
            0.95,
            chr(ord("A") + panel_index),
            transform=axis.transAxes,
            fontsize=SUBPLOT_LABEL_FONT_SIZE,
            fontweight="bold",
            ha="right",
            va="top",
            clip_on=True,
            zorder=10,
        )
        axis.axhline(
            0,
            color="black",
            linestyle="--",
            linewidth=1.0,
            zorder=1,
        )

        for peak_name in TARGET_MASSES:
            peak_results = metric_results.loc[
                metric_results["Peak"] == peak_name
            ].sort_values("Injection Number")
            if peak_results.empty:
                continue
            color = peak_colors[peak_name]
            plotted = axis.errorbar(
                peak_results["Injection Number"],
                peak_results["Delta"],
                yerr=peak_results["SEM"],
                marker="o",
                linestyle="-",
                linewidth=0.7,
                markersize=1.8,
                capsize=1.5,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color,
                elinewidth=0.6,
                capthick=0.6,
                label=peak_name,
                zorder=2,
            )
            if panel_index == 0:
                legend_handles.append(plotted)
                legend_labels.append(peak_name)

        set_delta_axis_limits(axis, metric_results, metric["minimum_limit"])
        axis.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        axis.set_xlabel("Number of Injections")
        axis.set_ylabel(metric["ylabel"])
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_linewidth(AXIS_STROKE_WIDTH)
        axis.spines["bottom"].set_linewidth(AXIS_STROKE_WIDTH)
        axis.tick_params(width=AXIS_STROKE_WIDTH)

    if legend_handles:
        figure.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=5,
            frameon=False,
            fontsize=LEGEND_FONT_SIZE,
            columnspacing=1.2,
            handlelength=1.8,
        )

    figure.tight_layout(rect=(0, 0, 1, 0.80), w_pad=1.5)
    return figure


def main() -> None:
    """Calculate Rabbit IgG deltas and render Figure S14."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": BASE_FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    measurements = collect_measurements()
    results = create_delta_results(measurements)
    results.to_csv(OUTPUT_CSV, index=False, float_format="%.10g")

    figure = make_figure(results)
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)

    maximum_measurements = (
        len(REPLICATE_FOLDERS) * len(EXPECTED_INJECTIONS) * len(TARGET_MASSES)
    )
    print(
        f"Matched {len(measurements)} of {maximum_measurements} possible "
        "Rabbit IgG peak measurements"
    )
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")
    print(f"Saved {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
