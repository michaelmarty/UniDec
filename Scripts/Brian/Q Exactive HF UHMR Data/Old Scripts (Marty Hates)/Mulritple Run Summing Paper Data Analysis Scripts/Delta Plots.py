"""
Plot changes in UniDec ion statistics relative to injection 21.

Mass centroid, FWHM, and sum-normalized peak area fraction are shown in three
panels. Every requested peak is overlaid in each panel with a consistent color.
Points are differences between replicate means and the injection 21 mean:

    delta = injection mean - injection 21 mean

Error bars are SEM in the original units of each metric.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator


DEFAULT_MASS_TOLERANCE = 50.0
DEFAULT_OUTPUT_FILENAME = "delta_plots.pdf"
REFERENCE_INJECTION = 21
MASS_COLUMN_INDEX = 0
AREA_COLUMN_INDEX = 5
MASS_CENTROID_COLUMN_INDEX = 6
FWHM_COLUMN_INDEX = 7
BASE_FONT_SIZE = 16
PANEL_TITLE_FONT_SIZE = 20
SUBPLOT_LABEL_FONT_SIZE = 19
LEGEND_FONT_SIZE = 13

METRICS = {
    "centroid": {
        "title": "Mass Centroid",
        "table_name": "Mass Centroid",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION} (Da)",
    },
    "fwhm": {
        "title": "FWHM",
        "table_name": "FWHM",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION} (Da)",
    },
    "area_fraction": {
        "title": "Sum-Normalized\nPeak Area Fraction",
        "table_name": "Sum-Normalized Peak Area Fraction",
        "ylabel": f"Delta from Injection {REFERENCE_INJECTION}",
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
        if output_filename.parent != Path(""):
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
    print("DELTA PLOTS")
    print("=" * 80)
    print(
        f"All values will be plotted as deltas from injection {REFERENCE_INJECTION}."
    )

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

    print("\n" + "=" * 80)
    print("CONFIGURATION SUMMARY")
    print("=" * 80)
    print(f"Replicate folders: {len(folders)}")
    for index, folder in enumerate(folders, 1):
        print(f"  {index}. {folder}")
    print(f"Output file: {Path(output_folder) / output_filename}")
    print(f"Reference injection: {REFERENCE_INJECTION}")
    print(f"Mass tolerance: {mass_tolerance:g} Da")
    print("Metrics: Mass centroid, FWHM, and sum-normalized peak area fraction")
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
    }


def read_dat_file(filepath):
    """Read one UniDec peak-parameter DAT file."""
    try:
        return pd.read_csv(filepath, sep=r"\s+", encoding="utf-8")
    except Exception as error:
        print(f"Error reading {filepath}: {error}")
        return None


def get_peak_measurements(df, filepath):
    """Return numeric mass, area, centroid, and FWHM columns from a DAT file."""
    required_columns = max(
        MASS_COLUMN_INDEX,
        AREA_COLUMN_INDEX,
        MASS_CENTROID_COLUMN_INDEX,
        FWHM_COLUMN_INDEX,
    ) + 1
    if df.shape[1] < required_columns:
        print(
            f"Warning: Expected at least {required_columns} columns in {filepath}; "
            f"found {df.shape[1]}"
        )
        return pd.DataFrame(columns=["mass", "area", "centroid", "fwhm"])

    measurements = pd.DataFrame(
        {
            "mass": pd.to_numeric(df.iloc[:, MASS_COLUMN_INDEX], errors="coerce"),
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


def match_mass_to_target(mass, target_masses, tolerance):
    """Match the DAT file's first-column mass to the closest target."""
    candidates = [
        (abs(mass - target_mass), name, target_mass)
        for name, target_mass in target_masses.items()
        if abs(mass - target_mass) <= tolerance
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

    Peak matching and tolerance use the first DAT column. If multiple measured
    peaks match a target, the closest one is retained. Area is normalized by the
    sum of all retained target-peak areas in that DAT file.
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
                name, target_mass = match_mass_to_target(
                    float(peak.mass), target_masses, mass_tolerance
                )
                if name is None:
                    continue
                difference = abs(float(peak.mass) - target_mass)
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
    """Calculate the raw mean and SEM for each peak, metric, and injection."""
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


def get_injection_values(data, injection_number):
    """Return the mean and SEM for a requested injection, if present."""
    for injection, mean, error in zip(
        data["injections"], data["mean"], data["se"]
    ):
        if injection == injection_number:
            return float(mean), float(error)
    return None, None


def calculate_deltas(
    raw_results, target_masses, reference_injection=REFERENCE_INJECTION
):
    """
    Subtract the reference-injection mean from every injection mean.

    SEM remains in the original units of the metric.
    """
    delta_results = {
        peak_name: {
            metric: {
                "injections": [],
                "mean": [],
                "se": [],
                "reference_mean": None,
            }
            for metric in METRICS
        }
        for peak_name in target_masses
    }

    for peak_name in target_masses:
        for metric in METRICS:
            raw_data = raw_results[peak_name][metric]
            reference_mean, _ = get_injection_values(
                raw_data, reference_injection
            )
            if reference_mean is None:
                print(
                    f"Warning: No injection {reference_injection} reference for "
                    f"{peak_name} ({METRICS[metric]['table_name']}); series omitted."
                )
                continue

            delta_data = delta_results[peak_name][metric]
            delta_data["reference_mean"] = reference_mean
            for injection, mean, error in zip(
                raw_data["injections"], raw_data["mean"], raw_data["se"]
            ):
                delta_data["injections"].append(injection)
                delta_data["mean"].append(float(mean) - reference_mean)
                delta_data["se"].append(float(error))

    return delta_results


def metric_has_data(results, target_masses, metric):
    """Return True if at least one peak has plottable data for a metric."""
    return any(results[peak_name][metric]["injections"] for peak_name in target_masses)


def set_delta_axis_limits(ax, results, target_masses, metric):
    """Set a symmetric delta axis that includes all SEM bars and zero."""
    bounds = [0.0]
    for peak_name in target_masses:
        data = results[peak_name][metric]
        if not data["mean"]:
            continue
        means = np.asarray(data["mean"], dtype=float)
        errors = np.asarray(data["se"], dtype=float)
        bounds.extend((means - errors).tolist())
        bounds.extend((means + errors).tolist())

    maximum = max(abs(value) for value in bounds if np.isfinite(value))
    if maximum == 0:
        maximum = 1.0
    limit = maximum * 1.10
    ax.set_ylim(-limit, limit)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, symmetric=True))


def get_peak_colors(target_masses):
    """Assign each peak a distinct, stable color."""
    peak_names = list(target_masses)
    if len(peak_names) <= 10:
        colors = plt.get_cmap("tab10").colors
    elif len(peak_names) <= 20:
        colors = plt.get_cmap("tab20").colors
    else:
        color_map = plt.get_cmap("gist_ncar")
        denominator = max(len(peak_names) - 1, 1)
        return {
            peak_name: color_map(index / denominator)
            for index, peak_name in enumerate(peak_names)
        }
    return {peak_name: colors[index] for index, peak_name in enumerate(peak_names)}


def plot_results(results, target_masses, output_folder, output_filename):
    """Create three panels, with all peaks overlaid in every metric panel."""
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = BASE_FONT_SIZE
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    metrics = list(METRICS)
    fig, axes = plt.subplots(
        1,
        len(metrics),
        squeeze=False,
        sharex=True,
        figsize=(17.5, 5.8),
    )
    axes = axes[0]
    peak_colors = get_peak_colors(target_masses)
    legend_handles = []
    legend_labels = []

    for panel_index, (ax, metric) in enumerate(zip(axes, metrics)):
        ax.set_title(
            METRICS[metric]["title"],
            fontsize=PANEL_TITLE_FONT_SIZE,
            fontweight="bold",
            pad=14,
        )
        ax.text(
            0.97,
            0.95,
            chr(ord("A") + panel_index),
            transform=ax.transAxes,
            fontsize=SUBPLOT_LABEL_FONT_SIZE,
            fontweight="bold",
            ha="right",
            va="top",
            clip_on=True,
            zorder=10,
        )
        ax.axhline(0, color="black", linestyle="--", linewidth=1.2, zorder=1)

        for peak_name in target_masses:
            data = results[peak_name][metric]
            if not data["injections"]:
                continue
            color = peak_colors[peak_name]
            plotted = ax.errorbar(
                data["injections"],
                data["mean"],
                yerr=data["se"],
                marker="o",
                linestyle="-",
                linewidth=1.5,
                markersize=4.5,
                capsize=3,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color,
                elinewidth=1.2,
                label=peak_name,
                zorder=2,
            )
            if panel_index == 0:
                legend_handles.append(plotted)
                legend_labels.append(peak_name)

        set_delta_axis_limits(ax, results, target_masses, metric)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        ax.set_xlabel("Number of Injections", fontweight="bold")
        ax.set_ylabel(METRICS[metric]["ylabel"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if not metric_has_data(results, target_masses, metric):
            ax.text(
                0.5,
                0.5,
                "No reference data available",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=min(len(legend_labels), 6),
            frameon=False,
            fontsize=LEGEND_FONT_SIZE,
        )

    fig.tight_layout(rect=(0, 0, 1, 0.91 if legend_handles else 1), w_pad=2.0)

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
        "Format: PDF (300 DPI, transparent, 17.5 x 5.8 inches, editable text)"
    )


def create_summary_tables(results, target_masses):
    """Print delta and SEM result tables."""
    print("\n" + "=" * 80)
    print(f"DELTAS FROM INJECTION {REFERENCE_INJECTION}")
    print("=" * 80)
    for peak_name in target_masses:
        for metric in METRICS:
            data = results[peak_name][metric]
            print(f"\n{peak_name} - {METRICS[metric]['table_name']}")
            if data["reference_mean"] is None:
                print(f"No injection {REFERENCE_INJECTION} reference available.")
                continue
            print(f"Injection {REFERENCE_INJECTION} mean: {data['reference_mean']:.6g}")
            print(f"{'Injections':<12} {'Delta':<18} {'SEM':<18}")
            print("-" * 48)
            for injection, mean, error in zip(
                data["injections"], data["mean"], data["se"]
            ):
                print(f"{injection:<12} {mean:<18.6f} {error:<18.6f}")


def save_summary_csv(results, target_masses, output_folder, output_filename):
    """Save plotted deltas and SEM values as one CSV."""
    rows = []
    for peak_name in target_masses:
        for metric in METRICS:
            data = results[peak_name][metric]
            for injection, mean, error in zip(
                data["injections"], data["mean"], data["se"]
            ):
                rows.append(
                    {
                        "Peak": peak_name,
                        "Metric": METRICS[metric]["table_name"],
                        "Reference Injection": REFERENCE_INJECTION,
                        "Reference Mean": data["reference_mean"],
                        "Injections": injection,
                        "Delta": mean,
                        "SEM": error,
                    }
                )

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    csv_filename = f"{Path(output_filename).stem}_summary.csv"
    csv_file = output_path / csv_filename
    columns = [
        "Peak",
        "Metric",
        "Reference Injection",
        "Reference Mean",
        "Injections",
        "Delta",
        "SEM",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(csv_file, index=False)
    print(f"\nSummary CSV saved as '{csv_file}'")


def main():
    """Run the interactive delta-plot workflow."""
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
    if REFERENCE_INJECTION not in all_data:
        print(
            f"\nInjection {REFERENCE_INJECTION} was not found. "
            "Deltas cannot be calculated."
        )
        return

    raw_results = calculate_statistics(
        all_data, settings["target_masses"], len(settings["folders"])
    )
    delta_results = calculate_deltas(
        raw_results, settings["target_masses"]
    )
    if not any(
        metric_has_data(delta_results, settings["target_masses"], metric)
        for metric in METRICS
    ):
        print(
            f"\nNo peaks had usable injection {REFERENCE_INJECTION} reference data."
        )
        return

    create_summary_tables(delta_results, settings["target_masses"])
    save_summary_csv(
        delta_results,
        settings["target_masses"],
        settings["output_folder"],
        settings["output_filename"],
    )
    plot_results(
        delta_results,
        settings["target_masses"],
        settings["output_folder"],
        settings["output_filename"],
    )
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
