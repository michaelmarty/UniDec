"""
Plot UniDec ion statistics as a peak-by-metric grid.

Rows are user-specified peaks. Columns are mass centroid, FWHM, and,
optionally, sum-normalized peak area fraction. Values are means +/- SEM
across replicate folders for each injection number.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, StrMethodFormatter


DEFAULT_MASS_TOLERANCE = 50.0
DEFAULT_OUTPUT_FILENAME = "ion_statistics_grid.pdf"
HIGHLIGHT_INJECTION = 21
AREA_COLUMN_INDEX = 5
MASS_CENTROID_COLUMN_INDEX = 6
FWHM_COLUMN_INDEX = 7
Y_TICK_INTERVALS_PER_SIDE = 2
BASE_FONT_SIZE = 16
GRID_HEADER_FONT_SIZE = 20
SUBPLOT_LABEL_FONT_SIZE = 19
LEGEND_FONT_SIZE = 13
FIGURE_WIDTH_PER_COLUMN = 5.0
FIGURE_HEIGHT_PER_ROW = 3.4

METRICS = {
    "centroid": {
        "title": "Mass Centroid",
        "ylabel": "Mass Centroid (Da)",
        "format": "{x:.0f}",
    },
    "fwhm": {
        "title": "FWHM",
        "ylabel": "FWHM (Da)",
        "format": "{x:.0f}",
    },
    "area_fraction": {
        "title": "Sum-Normalized\nPeak Area Fraction",
        "ylabel": "Area Fraction",
        "format": "{x:.2f}",
    },
}


def clean_path(path_string):
    """Normalize a path pasted from Windows Explorer or a terminal."""
    path_string = path_string.strip()
    if (
        len(path_string) >= 2
        and path_string[0] == path_string[-1]
        and path_string[0] in "\"'"
    ):
        path_string = path_string[1:-1]
    return path_string.strip()


def prompt_for_positive_int(prompt):
    """Prompt until the user enters a positive integer."""
    while True:
        try:
            value = int(input(prompt))
            if value > 0:
                return value
        except ValueError:
            pass
        print("Please enter a positive integer.")


def prompt_for_positive_float(prompt, default=None):
    """Prompt until the user enters a positive number."""
    while True:
        raw_value = input(prompt).strip()
        if not raw_value and default is not None:
            return float(default)
        try:
            value = float(raw_value)
            if value > 0:
                return value
        except ValueError:
            pass
        print("Please enter a positive number.")


def prompt_for_yes_no(prompt, default=None):
    """Prompt for a yes/no response and return a bool."""
    suffix = ""
    if default is True:
        suffix = " [Y/n]"
    elif default is False:
        suffix = " [y/N]"

    while True:
        answer = input(f"{prompt}{suffix}: ").strip().lower()
        if not answer and default is not None:
            return default
        if answer in ("yes", "y"):
            return True
        if answer in ("no", "n"):
            return False
        print("Please enter 'yes' or 'no'.")


def prompt_for_paths(count, prompt_name):
    """Prompt for a fixed number of paths."""
    paths = []
    print(f"\nEnter the paths for {count} {prompt_name}:")
    print("(Quotes around pasted paths are handled automatically.)")
    for index in range(count):
        while True:
            path = clean_path(input(f"  {prompt_name[:-1].title()} {index + 1} path: "))
            if path:
                paths.append(path)
                break
            print("  Please enter a valid path.")
    return paths


def prompt_for_output_filename(default_filename):
    """Prompt for a PDF file name without a directory component."""
    while True:
        filename = clean_path(
            input(f"\nEnter output PDF file name [default: {default_filename}]: ")
        )
        if not filename:
            return default_filename

        output_filename = Path(filename)
        if output_filename.parent != Path("."):
            print("  Please enter a file name only, not a folder path.")
            continue

        filename = output_filename.name
        if not filename.lower().endswith(".pdf"):
            filename = f"{filename.rstrip('.')}.pdf"
        if filename != ".pdf":
            return filename
        print("  Please enter a valid file name.")


def get_user_inputs():
    """Collect and confirm all analysis settings."""
    print("=" * 80)
    print("ION STATISTICS GRID")
    print("=" * 80)

    num_replicates = prompt_for_positive_int("\nEnter the number of replicate folders: ")
    folders = prompt_for_paths(num_replicates, "replicate folders")

    print("\nEnter the output folder path (leave blank for the current directory):")
    output_folder = clean_path(input("  Output folder: ")) or "."
    output_filename = prompt_for_output_filename(DEFAULT_OUTPUT_FILENAME)

    num_peaks = prompt_for_positive_int("\nEnter the number of peaks to track: ")
    target_masses = {}
    print(f"\nEnter {num_peaks} expected masses and peak names:")
    for index in range(num_peaks):
        print(f"\n  Peak {index + 1}:")
        while True:
            name = input("    Peak name: ").strip()
            if not name:
                print("    Please enter a peak name.")
            elif name in target_masses:
                print("    Peak names must be unique.")
            else:
                break
        target_masses[name] = prompt_for_positive_float("    Expected mass (Da): ")

    mass_tolerance = prompt_for_positive_float(
        f"\nEnter mass tolerance for matching (Da) "
        f"[default: {DEFAULT_MASS_TOLERANCE:g}]: ",
        DEFAULT_MASS_TOLERANCE,
    )
    show_theoretical_mass = prompt_for_yes_no(
        "\nDisplay the theoretical mass line in the mass centroid plots", True
    )
    include_area_fraction = prompt_for_yes_no(
        "Include sum-normalized peak area fraction plots", True
    )

    print("\n" + "=" * 80)
    print("CONFIGURATION SUMMARY")
    print("=" * 80)
    print(f"Replicate folders: {len(folders)}")
    for index, folder in enumerate(folders, 1):
        print(f"  {index}. {folder}")
    print(f"Output file: {Path(output_folder) / output_filename}")
    print(f"Mass tolerance: {mass_tolerance:g} Da")
    print(f"Theoretical mass line: {'Yes' if show_theoretical_mass else 'No'}")
    print(
        "Sum-normalized peak area fraction: "
        f"{'Included' if include_area_fraction else 'Excluded'}"
    )
    print("Peaks:")
    for name, mass in target_masses.items():
        print(f"  {name}: {mass:.2f} Da")

    if not prompt_for_yes_no("\nProceed with this configuration"):
        print("\nRestarting configuration...")
        return get_user_inputs()

    return {
        "folders": folders,
        "output_folder": output_folder,
        "output_filename": output_filename,
        "target_masses": target_masses,
        "mass_tolerance": mass_tolerance,
        "show_theoretical_mass": show_theoretical_mass,
        "include_area_fraction": include_area_fraction,
    }


def read_dat_file(filepath):
    """Read one UniDec peak-parameter DAT file."""
    try:
        return pd.read_csv(filepath, sep=r"\s+", encoding="utf-8")
    except Exception as error:
        print(f"Error reading {filepath}: {error}")
        return None


def get_peak_measurements(df, filepath):
    """Return numeric area, centroid, and FWHM columns from a DAT file."""
    required_columns = max(
        AREA_COLUMN_INDEX, MASS_CENTROID_COLUMN_INDEX, FWHM_COLUMN_INDEX
    ) + 1
    if df.shape[1] < required_columns:
        print(
            f"Warning: Expected at least {required_columns} columns in {filepath}; "
            f"found {df.shape[1]}"
        )
        return pd.DataFrame(columns=["area", "centroid", "fwhm"])

    measurements = pd.DataFrame(
        {
            "area": pd.to_numeric(df.iloc[:, AREA_COLUMN_INDEX], errors="coerce"),
            "centroid": pd.to_numeric(
                df.iloc[:, MASS_CENTROID_COLUMN_INDEX], errors="coerce"
            ),
            "fwhm": pd.to_numeric(df.iloc[:, FWHM_COLUMN_INDEX], errors="coerce"),
        }
    ).dropna()
    if measurements.empty:
        print(f"Warning: No complete numeric peak measurements found in {filepath}")
    return measurements


def match_centroid_to_target(mass_centroid, target_masses, tolerance):
    """Return the closest target name and mass within tolerance."""
    candidates = [
        (abs(mass_centroid - target_mass), name, target_mass)
        for name, target_mass in target_masses.items()
        if abs(mass_centroid - target_mass) <= tolerance
    ]
    if not candidates:
        return None, None
    _, name, target_mass = min(candidates)
    return name, target_mass


def extract_injection_number(filename):
    """Extract an injection number from a name such as '*merged_21_*.dat'."""
    match = re.search(r"merged_(\d+)_", filename)
    return int(match.group(1)) if match else None


def process_all_data(folders, target_masses, mass_tolerance):
    """
    Return {injection: {replicate: {peak: metric values}}}.

    When more than one measured peak matches a target, the closest centroid is
    retained. Area is normalized by the sum of all retained target-peak areas
    in that DAT file.
    """
    all_data = {}

    for replicate, folder in enumerate(folders, 1):
        folder_path = Path(folder)
        if not folder_path.exists():
            print(f"Warning: Folder does not exist: {folder}")
            continue

        dat_files = sorted(folder_path.glob("*_peakparam.dat"))
        print(f"\nProcessing Rep {replicate}: Found {len(dat_files)} DAT files")

        for dat_file in dat_files:
            injection = extract_injection_number(dat_file.name)
            if injection is None:
                print(f"Warning: Could not extract injection number from {dat_file.name}")
                continue

            df = read_dat_file(dat_file)
            if df is None or df.empty:
                continue
            measurements = get_peak_measurements(df, dat_file)
            if measurements.empty:
                continue

            matched = {}
            match_differences = {}
            for peak in measurements.itertuples(index=False):
                name, target_mass = match_centroid_to_target(
                    float(peak.centroid), target_masses, mass_tolerance
                )
                if name is None:
                    continue
                difference = abs(float(peak.centroid) - target_mass)
                if difference < match_differences.get(name, float("inf")):
                    matched[name] = {
                        "centroid": float(peak.centroid),
                        "fwhm": float(peak.fwhm),
                        "area": float(peak.area),
                    }
                    match_differences[name] = difference

            if not matched:
                continue

            total_area = sum(values["area"] for values in matched.values())
            if total_area <= 0:
                print(f"Warning: Matched peak areas sum to zero or less in {dat_file}")

            rep_data = all_data.setdefault(injection, {}).setdefault(replicate, {})
            for name, values in matched.items():
                rep_data[name] = {
                    "centroid": values["centroid"],
                    "fwhm": values["fwhm"],
                    "area_fraction": (
                        values["area"] / total_area if total_area > 0 else np.nan
                    ),
                }

    return all_data


def standard_error(values):
    """Calculate SEM; a single value has an SEM of zero."""
    return float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0


def calculate_statistics(all_data, target_masses, num_replicates):
    """Calculate mean and SEM for every peak, metric, and injection."""
    results = {
        peak_name: {
            metric: {"injections": [], "mean": [], "se": []}
            for metric in METRICS
        }
        for peak_name in target_masses
    }

    for injection in sorted(all_data):
        for peak_name in target_masses:
            for metric in METRICS:
                values = []
                for replicate in range(1, num_replicates + 1):
                    value = (
                        all_data[injection]
                        .get(replicate, {})
                        .get(peak_name, {})
                        .get(metric)
                    )
                    if value is not None and np.isfinite(value):
                        values.append(value)
                if values:
                    data = results[peak_name][metric]
                    data["injections"].append(injection)
                    data["mean"].append(float(np.mean(values)))
                    data["se"].append(standard_error(values))
    return results


def get_nice_tick_step(raw_step):
    """Round a tick interval up to a readable 1, 2, 5, or 10 multiple."""
    if raw_step <= 0:
        return 1.0
    magnitude = 10 ** np.floor(np.log10(raw_step))
    normalized = raw_step / magnitude
    for candidate in (1, 2, 5, 10):
        if normalized <= candidate:
            return float(candidate * magnitude)
    return float(10 * magnitude)


def set_metric_axis_limits(ax, data, metric, theoretical_mass=None):
    """Scale a metric axis to its values and optional theoretical mass."""
    values = []
    if data["mean"]:
        means = np.asarray(data["mean"], dtype=float)
        errors = np.asarray(data["se"], dtype=float)
        values.extend((means - errors).tolist())
        values.extend((means + errors).tolist())
    if theoretical_mass is not None:
        values.append(float(theoretical_mass))
    if not values:
        return

    lower_value = min(values)
    upper_value = max(values)
    span = upper_value - lower_value
    minimum_padding = 0.01 if metric == "area_fraction" else 1
    padding = max(
        span * 0.08,
        max(abs(lower_value), abs(upper_value), 1) * 0.00005,
        minimum_padding,
    )
    lower_limit = lower_value - padding
    upper_limit = upper_value + padding
    tick_step = get_nice_tick_step(
        (upper_limit - lower_limit) / (Y_TICK_INTERVALS_PER_SIDE * 2)
    )
    tick_start = np.floor(lower_limit / tick_step) * tick_step
    tick_end = np.ceil(upper_limit / tick_step) * tick_step
    ticks = np.arange(tick_start, tick_end + tick_step * 0.5, tick_step)
    ax.set_ylim(ticks[0], ticks[-1])
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(StrMethodFormatter(METRICS[metric]["format"]))
    ax.yaxis.get_offset_text().set_visible(False)


def get_highlight_band(data, injection_number):
    """Return the mean and SEM at one injection number."""
    for injection, mean, error in zip(
        data["injections"], data["mean"], data["se"]
    ):
        if injection == injection_number:
            return float(mean), float(error)
    return None, None


def get_subplot_label(index):
    """Return spreadsheet-style subplot labels: A ... Z, AA, AB, etc."""
    label = ""
    while True:
        index, remainder = divmod(index, 26)
        label = chr(ord("A") + remainder) + label
        if index == 0:
            return label
        index -= 1


def plot_results(
    results,
    target_masses,
    output_folder,
    output_filename,
    show_theoretical_mass,
    include_area_fraction,
):
    """Create a grid whose rows are peaks and columns are metrics."""
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = BASE_FONT_SIZE
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    metrics = ["centroid", "fwhm"]
    if include_area_fraction:
        metrics.append("area_fraction")

    peak_names = list(target_masses)
    figure_size = (
        FIGURE_WIDTH_PER_COLUMN * len(metrics) + 1.5,
        FIGURE_HEIGHT_PER_ROW * len(peak_names) + 1.5,
    )
    fig, axes = plt.subplots(
        len(peak_names),
        len(metrics),
        squeeze=False,
        sharex="col",
        figsize=figure_size,
    )
    point_color = plt.get_cmap("tab10")(0)
    highlight_color = plt.get_cmap("tab10")(2)

    for row, peak_name in enumerate(peak_names):
        for column, metric in enumerate(metrics):
            ax = axes[row, column]
            data = results[peak_name][metric]

            if row == 0:
                ax.set_title(
                    METRICS[metric]["title"],
                    fontsize=GRID_HEADER_FONT_SIZE,
                    fontweight="bold",
                    pad=14,
                )

            if column == 0:
                ax.annotate(
                    peak_name,
                    xy=(0, 0.5),
                    xytext=(-ax.yaxis.labelpad - 38, 0),
                    xycoords=ax.yaxis.label,
                    textcoords="offset points",
                    ha="right",
                    va="center",
                    fontsize=GRID_HEADER_FONT_SIZE,
                    fontweight="bold",
                    rotation=90,
                )

            subplot_index = row * len(metrics) + column
            ax.text(
                -0.14,
                1.05,
                get_subplot_label(subplot_index),
                transform=ax.transAxes,
                fontsize=SUBPLOT_LABEL_FONT_SIZE,
                fontweight="bold",
                ha="left",
                va="bottom",
                clip_on=False,
            )

            theoretical_for_limits = None
            if metric == "centroid" and show_theoretical_mass:
                theoretical_for_limits = target_masses[peak_name]
                ax.axhline(
                    theoretical_for_limits,
                    color="black",
                    linestyle="-",
                    linewidth=1.4,
                    label="Theoretical mass",
                    zorder=1.5,
                )

            highlight_mean, highlight_error = get_highlight_band(
                data, HIGHLIGHT_INJECTION
            )
            if highlight_mean is not None and highlight_error > 0:
                ax.axhspan(
                    highlight_mean - highlight_error,
                    highlight_mean + highlight_error,
                    facecolor=highlight_color,
                    alpha=0.45,
                    linewidth=0,
                    label=f"{HIGHLIGHT_INJECTION} injections mean +/- SEM",
                    zorder=1,
                )

            if data["injections"]:
                ax.errorbar(
                    data["injections"],
                    data["mean"],
                    yerr=data["se"],
                    marker="o",
                    linestyle="-",
                    linewidth=1.5,
                    markersize=4,
                    capsize=3,
                    color=point_color,
                    markerfacecolor=point_color,
                    markeredgecolor=point_color,
                    ecolor="black",
                    elinewidth=1.2,
                    zorder=2,
                )

            set_metric_axis_limits(ax, data, metric, theoretical_for_limits)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            ax.set_ylabel(METRICS[metric]["ylabel"])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row == len(peak_names) - 1:
                ax.set_xlabel("Number of Injections", fontweight="bold")

            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(frameon=False, fontsize=LEGEND_FONT_SIZE, loc="best")

    fig.tight_layout(rect=(0.03, 0, 1, 1), h_pad=2.2, w_pad=2.0)

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / output_filename
    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        format="pdf",
    )
    plt.show()

    print(f"\nPlot saved as '{output_file}'")
    print(
        "Format: PDF "
        f"(300 DPI, transparent, {figure_size[0]:g} x {figure_size[1]:g} inches, "
        "editable text)"
    )


def create_summary_tables(results, target_masses, include_area_fraction):
    """Print compact result tables for the plotted metrics."""
    metrics = ["centroid", "fwhm"]
    if include_area_fraction:
        metrics.append("area_fraction")

    print("\n" + "=" * 80)
    print("SUMMARY TABLES")
    print("=" * 80)
    for peak_name in target_masses:
        for metric in metrics:
            print(f"\n{peak_name} - {METRICS[metric]['title']}")
            print(f"{'Injections':<12} {'Mean':<18} {'SEM':<18}")
            print("-" * 48)
            data = results[peak_name][metric]
            for injection, mean, error in zip(
                data["injections"], data["mean"], data["se"]
            ):
                print(f"{injection:<12} {mean:<18.4f} {error:<18.4f}")


def save_summary_csv(
    results,
    target_masses,
    include_area_fraction,
    output_folder,
    output_filename,
):
    """Save the printed injection, mean, and SEM results as one CSV."""
    metrics = ["centroid", "fwhm"]
    if include_area_fraction:
        metrics.append("area_fraction")

    rows = []
    for peak_name in target_masses:
        for metric in metrics:
            data = results[peak_name][metric]
            metric_name = METRICS[metric]["title"].replace("\n", " ")
            for injection, mean, error in zip(
                data["injections"], data["mean"], data["se"]
            ):
                rows.append(
                    {
                        "Peak": peak_name,
                        "Metric": metric_name,
                        "Injections": injection,
                        "Mean": mean,
                        "SEM": error,
                    }
                )

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    csv_filename = f"{Path(output_filename).stem}_summary.csv"
    csv_file = output_path / csv_filename
    pd.DataFrame(
        rows,
        columns=["Peak", "Metric", "Injections", "Mean", "SEM"],
    ).to_csv(csv_file, index=False)
    print(f"\nSummary CSV saved as '{csv_file}'")


def main():
    """Run the interactive ion-statistics workflow."""
    settings = get_user_inputs()
    print("\nStarting DAT file analysis...")
    all_data = process_all_data(
        settings["folders"],
        settings["target_masses"],
        settings["mass_tolerance"],
    )
    if not all_data:
        print("\nNo data was processed. Please check folder paths and file names.")
        return

    results = calculate_statistics(
        all_data, settings["target_masses"], len(settings["folders"])
    )
    create_summary_tables(
        results,
        settings["target_masses"],
        settings["include_area_fraction"],
    )
    save_summary_csv(
        results,
        settings["target_masses"],
        settings["include_area_fraction"],
        settings["output_folder"],
        settings["output_filename"],
    )
    plot_results(
        results,
        settings["target_masses"],
        settings["output_folder"],
        settings["output_filename"],
        settings["show_theoretical_mass"],
        settings["include_area_fraction"],
    )
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
