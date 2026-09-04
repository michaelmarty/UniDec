"""Create the IgM-like fusion ion-statistics grid used in Figure S16."""

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
import numpy as np
import pandas as pd


REPLICATE_FOLDERS = [
    Path(r"Z:\Group Share\BHT\Temp\Paper\Stats\Amgen IgM-like Fusion Rep 1"),
    Path(r"Z:\Group Share\BHT\Temp\Paper\Stats\Amgen IgM-like Fusion Rep 2"),
    Path(r"Z:\Group Share\BHT\Temp\Paper\Stats\Amgen IgM-like Fusion Rep 3"),
]

TARGET_MASSES = {
    "IgM-like Fusion Protein": 1183000.0,
}

MASS_TOLERANCE_DA = 100000.0
EXPECTED_INJECTIONS = tuple(range(1, 22))
HIGHLIGHT_INJECTION = 21
OUTPUT_DPI = 300
OUTPUT_PADDING_INCHES = 0.03
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S16"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S16.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S16.pdf"
OUTPUT_CSV = OUTPUT_DIRECTORY / "Figure S16_results.csv"

Y_TICK_INTERVALS_PER_SIDE = 2
BASE_FONT_SIZE = 10
GRID_HEADER_FONT_SIZE = 10
SUBPLOT_LABEL_FONT_SIZE = 11
LEGEND_FONT_SIZE = 10
LEGEND_LOCATION = "lower left"
FIGURE_SIZE_INCHES = (6.5, 3.0)
PEAK_LABEL_OFFSET_POINTS = 10
AXIS_STROKE_WIDTH = 0.8
MEASUREMENT_LABEL = "fusion-protein peak"

REQUIRED_COLUMNS = {
    "Mass",
    "MassCentroid",
    "MassFWHM",
}

METRICS = {
    "centroid": {
        "title": "Mass Centroid",
        "ylabel": "Mass Centroid (Da)",
        "source": "Mass Centroid (Da)",
        "mean": "Mean Mass Centroid (Da)",
        "sem": "SEM Mass Centroid (Da)",
        "replicate": "Replicate {replicate} Mass Centroid (Da)",
        "format": "{x:,.0f}",
    },
    "fwhm": {
        "title": "FWHM",
        "ylabel": "FWHM (Da)",
        "source": "FWHM (Da)",
        "mean": "Mean FWHM (Da)",
        "sem": "SEM FWHM (Da)",
        "replicate": "Replicate {replicate} FWHM (Da)",
        "format": "{x:,.0f}",
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


def match_target_peak(
    table: pd.DataFrame,
    target_mass: float,
    path: Path,
) -> pd.Series | None:
    """Return the closest peak whose mass centroid is within tolerance."""
    mass_error = (table["MassCentroid"] - target_mass).abs()
    candidates = table.loc[mass_error <= MASS_TOLERANCE_DA]
    if candidates.empty:
        return None
    closest_index = mass_error.loc[candidates.index].idxmin()
    return table.loc[closest_index]


def collect_measurements() -> pd.DataFrame:
    """Collect the matched IgM-like fusion peak from all replicate DAT files."""
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
            matched_peaks = {
                protein: match_target_peak(table, target_mass, path)
                for protein, target_mass in TARGET_MASSES.items()
            }
            available_peaks = {
                protein: peak
                for protein, peak in matched_peaks.items()
                if peak is not None
            }

            for protein, peak in matched_peaks.items():
                if peak is None:
                    print(
                        f"Warning: No {protein} peak within "
                        f"{MASS_TOLERANCE_DA:g} Da in {path.name}"
                    )
            if not available_peaks:
                continue

            for protein, target_mass in TARGET_MASSES.items():
                peak = matched_peaks[protein]
                if peak is None:
                    continue
                records.append(
                    {
                        "Protein": protein,
                        "Theoretical Mass (Da)": target_mass,
                        "Injection Number": injection_number,
                        "Replicate": replicate,
                        "Matched Mass (Da)": float(peak["Mass"]),
                        "Mass Centroid (Da)": float(peak["MassCentroid"]),
                        "FWHM (Da)": float(peak["MassFWHM"]),
                        "Source File": str(path),
                    }
                )

    measurements = pd.DataFrame.from_records(records)
    maximum_rows = (
        len(REPLICATE_FOLDERS) * len(EXPECTED_INJECTIONS) * len(TARGET_MASSES)
    )
    if measurements.empty or len(measurements) > maximum_rows:
        raise ValueError(
            f"Expected between 1 and {maximum_rows} matched measurements; "
            f"found {len(measurements)}"
        )
    return measurements


def standard_error(values: np.ndarray) -> float:
    """Return the standard error of the mean for replicate values."""
    if values.size < 2:
        return 0.0
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def create_results_table(measurements: pd.DataFrame) -> pd.DataFrame:
    """Create the wide replicate, mean, and SEM table used for plotting."""
    rows: list[dict[str, float | int | str]] = []

    for protein, theoretical_mass in TARGET_MASSES.items():
        protein_data = measurements.loc[measurements["Protein"] == protein]
        for injection_number in EXPECTED_INJECTIONS:
            group = protein_data.loc[
                protein_data["Injection Number"] == injection_number
            ].sort_values("Replicate")
            if group.empty:
                raise ValueError(
                    f"No replicate data for {protein}, injection {injection_number}"
                )

            row: dict[str, float | int | str] = {
                "Protein": protein,
                "Theoretical Mass (Da)": theoretical_mass,
                "Injection Number": injection_number,
            }
            for metric in METRICS.values():
                values = group[metric["source"]].to_numpy(dtype=float)
                for replicate in range(1, len(REPLICATE_FOLDERS) + 1):
                    replicate_data = group.loc[group["Replicate"] == replicate]
                    value = (
                        float(replicate_data.iloc[0][metric["source"]])
                        if not replicate_data.empty
                        else np.nan
                    )
                    row[metric["replicate"].format(replicate=replicate)] = value
                row[metric["mean"]] = float(np.mean(values))
                row[metric["sem"]] = standard_error(values)
            rows.append(row)

    return pd.DataFrame.from_records(rows)


def get_nice_tick_step(raw_step: float) -> float:
    """Round a tick interval up to a readable 1, 2, 5, or 10 multiple."""
    if raw_step <= 0:
        return 1.0
    magnitude = 10 ** np.floor(np.log10(raw_step))
    normalized = raw_step / magnitude
    for candidate in (1, 2, 5, 10):
        if normalized <= candidate:
            return float(candidate * magnitude)
    return float(10 * magnitude)


def set_metric_axis_limits(
    axis: plt.Axes,
    means: np.ndarray,
    errors: np.ndarray,
    metric_name: str,
) -> None:
    """Scale an axis using the nice-tick logic from Ion Statistics Grid.py."""
    values = np.concatenate((means - errors, means + errors))

    lower_value = float(np.min(values))
    upper_value = float(np.max(values))
    span = upper_value - lower_value
    padding = max(
        span * 0.08,
        max(abs(lower_value), abs(upper_value), 1.0) * 0.00005,
        1.0,
    )
    lower_limit = lower_value - padding
    upper_limit = upper_value + padding
    tick_step = get_nice_tick_step(
        (upper_limit - lower_limit) / (Y_TICK_INTERVALS_PER_SIDE * 2)
    )
    tick_start = np.floor(lower_limit / tick_step) * tick_step
    tick_end = np.ceil(upper_limit / tick_step) * tick_step
    ticks = np.arange(tick_start, tick_end + tick_step * 0.5, tick_step)
    axis.set_ylim(ticks[0], ticks[-1])
    axis.set_yticks(ticks)
    axis.yaxis.set_major_formatter(StrMethodFormatter(METRICS[metric_name]["format"]))
    axis.yaxis.get_offset_text().set_visible(False)


def subplot_label(index: int) -> str:
    """Return spreadsheet-style subplot labels: A ... Z, AA, AB, etc."""
    label = ""
    while True:
        index, remainder = divmod(index, 26)
        label = chr(ord("A") + remainder) + label
        if index == 0:
            return label
        index -= 1


def make_figure(results: pd.DataFrame) -> plt.Figure:
    """Build the grid using the plotting logic from Ion Statistics Grid.py."""
    figure, axes = plt.subplots(
        nrows=len(TARGET_MASSES),
        ncols=len(METRICS),
        squeeze=False,
        sharex="col",
        figsize=FIGURE_SIZE_INCHES,
    )
    point_color = plt.get_cmap("tab10")(0)
    highlight_color = plt.get_cmap("tab10")(2)

    for row, (protein, theoretical_mass) in enumerate(TARGET_MASSES.items()):
        protein_results = results.loc[results["Protein"] == protein].sort_values(
            "Injection Number"
        )
        injection_numbers = protein_results["Injection Number"].to_numpy(dtype=float)

        for column, (metric_name, metric) in enumerate(METRICS.items()):
            axis = axes[row, column]
            means = protein_results[metric["mean"]].to_numpy(dtype=float)
            errors = protein_results[metric["sem"]].to_numpy(dtype=float)

            if row == 0:
                axis.set_title(
                    metric["title"],
                    fontsize=GRID_HEADER_FONT_SIZE,
                    fontweight="bold",
                    pad=14,
                )

            if column == 0:
                axis.annotate(
                    protein,
                    xy=(0, 0.5),
                    xytext=(-axis.yaxis.labelpad - PEAK_LABEL_OFFSET_POINTS, 0),
                    xycoords=axis.yaxis.label,
                    textcoords="offset points",
                    ha="right",
                    va="center",
                    fontsize=GRID_HEADER_FONT_SIZE,
                    fontweight="bold",
                    rotation=90,
                )

            axis.text(
                0.97,
                0.95,
                subplot_label(row * len(METRICS) + column),
                transform=axis.transAxes,
                fontsize=SUBPLOT_LABEL_FONT_SIZE,
                fontweight="bold",
                ha="right",
                va="top",
                clip_on=True,
                zorder=10,
            )

            highlight_row = protein_results.loc[
                protein_results["Injection Number"] == HIGHLIGHT_INJECTION
            ]
            if not highlight_row.empty:
                highlight_mean = float(highlight_row.iloc[0][metric["mean"]])
                highlight_error = float(highlight_row.iloc[0][metric["sem"]])
                if highlight_error > 0:
                    axis.axhspan(
                        highlight_mean - highlight_error,
                        highlight_mean + highlight_error,
                        facecolor=highlight_color,
                        alpha=0.45,
                        linewidth=0,
                        label=f"{HIGHLIGHT_INJECTION} INJ ± SEM",
                        zorder=1,
                    )

            axis.errorbar(
                injection_numbers,
                means,
                yerr=errors,
                marker="o",
                linestyle="-",
                linewidth=1.5,
                markersize=2.5,
                capsize=3,
                color=point_color,
                markerfacecolor=point_color,
                markeredgecolor=point_color,
                ecolor="black",
                elinewidth=AXIS_STROKE_WIDTH,
                capthick=AXIS_STROKE_WIDTH,
                zorder=2,
            )

            set_metric_axis_limits(
                axis,
                means,
                errors,
                metric_name,
            )
            axis.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            axis.set_ylabel(metric["ylabel"])
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.spines["left"].set_linewidth(AXIS_STROKE_WIDTH)
            axis.spines["bottom"].set_linewidth(AXIS_STROKE_WIDTH)
            axis.tick_params(width=AXIS_STROKE_WIDTH)

            if row == len(TARGET_MASSES) - 1:
                axis.set_xlabel("Number of injections")

            handles, _ = axis.get_legend_handles_labels()
            if handles and row == 0 and column == 0:
                axis.legend(
                    frameon=False,
                    fontsize=LEGEND_FONT_SIZE,
                    loc=LEGEND_LOCATION,
                )

    figure.tight_layout(rect=(0.03, 0, 1, 1), h_pad=2.2, w_pad=2.0)
    return figure


def main() -> None:
    """Calculate statistics and render Figure S16 plus its results CSV."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": BASE_FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    measurements = collect_measurements()
    results = create_results_table(measurements)
    results.to_csv(OUTPUT_CSV, index=False, float_format="%.10g")

    figure = make_figure(results)
    figure.savefig(
        OUTPUT_PNG,
        dpi=OUTPUT_DPI,
        transparent=True,
        bbox_inches="tight",
        pad_inches=OUTPUT_PADDING_INCHES,
    )
    figure.savefig(
        OUTPUT_PDF,
        dpi=OUTPUT_DPI,
        transparent=True,
        bbox_inches="tight",
        pad_inches=OUTPUT_PADDING_INCHES,
    )
    plt.close(figure)

    print(f"Matched {len(measurements)} {MEASUREMENT_LABEL} measurements")
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")
    print(f"Saved {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
