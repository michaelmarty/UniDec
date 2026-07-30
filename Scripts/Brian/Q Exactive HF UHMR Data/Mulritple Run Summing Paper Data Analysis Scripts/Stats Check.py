"""
Run one-way ANOVA and Dunnett post-hoc comparisons on UniDec *_peakparam.dat
files for mass centroid, normalized area ratio, normalized height ratio, and FWHM.
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


DEFAULT_MASS_TOLERANCE = 50.0
DEFAULT_CONTROL_INJECTION = 21
DEFAULT_ALPHA = 0.05

MEASUREMENT_COLUMNS = {
    "Height": 4,
    "Area": 5,
    "Mass Centroid": 6,
    "FWHM": 7,
}
STAT_METRICS = ["Mass Centroid", "Area Ratio", "Height Ratio", "FWHM"]

MASS_CENTROID_COLUMN_INDEX = MEASUREMENT_COLUMNS["Mass Centroid"]
REQUIRED_COLUMN_COUNT = max(MEASUREMENT_COLUMNS.values()) + 1


def clean_path(path_string):
    """
    Normalize a path pasted from Windows Explorer or the terminal.
    """
    path_string = path_string.strip()

    if (path_string.startswith('"') and path_string.endswith('"')) or \
            (path_string.startswith("'") and path_string.endswith("'")):
        path_string = path_string[1:-1]

    return path_string.strip()


def prompt_for_positive_int(prompt, default=None):
    """
    Prompt until the user enters a positive integer.
    """
    while True:
        raw_value = input(prompt).strip()
        if not raw_value and default is not None:
            return int(default)

        try:
            value = int(raw_value)
            if value > 0:
                return value
            print("Please enter a positive number.")
        except ValueError:
            print("Please enter a valid integer.")


def prompt_for_positive_float(prompt, default=None):
    """
    Prompt until the user enters a positive float.
    """
    while True:
        raw_value = input(prompt).strip()
        if not raw_value and default is not None:
            return float(default)

        try:
            value = float(raw_value)
            if value > 0:
                return value
            print("Please enter a positive value.")
        except ValueError:
            print("Please enter a valid number.")


def prompt_for_paths(count, prompt_name):
    """
    Prompt for a fixed number of folder paths.
    """
    paths = []
    print(f"\nEnter the paths for {count} {prompt_name}:")
    print("(Tip: Copy-pasted Windows paths and quotes will be handled automatically.)")

    for i in range(count):
        while True:
            path = clean_path(input(f"  {prompt_name[:-1].title()} {i + 1} path: "))
            if path:
                paths.append(path)
                break
            print("  Please enter a valid path.")

    return paths


def get_user_inputs():
    """
    Get all user inputs needed for the statistical analysis.
    """
    print("=" * 80)
    print("DAT FILE DUNNETT ANOVA STATISTICS")
    print("=" * 80)

    num_folders = prompt_for_positive_int("\nEnter the number of DAT file folders: ")
    folders = prompt_for_paths(num_folders, "DAT file folders")

    print("\nEnter the output folder path:")
    print("(Leave blank to use the current directory.)")
    output_folder = clean_path(input("  Output folder: "))
    if not output_folder:
        output_folder = "."
        print(f"  Using current directory: {output_folder}")

    num_masses = prompt_for_positive_int("\nEnter the number of peaks to track: ")
    target_masses = {}
    print(f"\nEnter {num_masses} expected mass values with their peak names:")

    for i in range(num_masses):
        print(f"\n  Mass {i + 1}:")
        while True:
            name = input("    Name (e.g., 'Cytochrome C'): ").strip()
            if name:
                break
            print("    Please enter a name.")

        target_masses[name] = prompt_for_positive_float("    Expected mass (Da): ")

    mass_tolerance = prompt_for_positive_float(
        f"\nEnter mass tolerance for matching (Da) [default: {DEFAULT_MASS_TOLERANCE:g}]: ",
        DEFAULT_MASS_TOLERANCE,
    )
    control_injection = prompt_for_positive_int(
        f"\nEnter the control injection number for Dunnett comparisons "
        f"[default: {DEFAULT_CONTROL_INJECTION}]: ",
        DEFAULT_CONTROL_INJECTION,
    )
    alpha = prompt_for_positive_float(
        f"\nEnter alpha/significance cutoff [default: {DEFAULT_ALPHA:g}]: ",
        DEFAULT_ALPHA,
    )

    print("\n" + "=" * 80)
    print("CONFIGURATION SUMMARY")
    print("=" * 80)
    print(f"\nNumber of DAT file folders: {num_folders}")
    print("\nDAT file folders:")
    for i, folder in enumerate(folders, 1):
        print(f"  {i}. {folder}")
    print(f"\nOutput folder: {output_folder}")
    print(f"Mass tolerance: {mass_tolerance} Da")
    print(f"Control injection: {control_injection}")
    print(f"Alpha: {alpha}")
    print("\nPeak targets:")
    for name, mass in target_masses.items():
        print(f"  {name}: {mass:.2f} Da")

    while True:
        confirm = input("\nProceed with this configuration? (yes/no): ").strip().lower()
        if confirm in ["yes", "y"]:
            break
        if confirm in ["no", "n"]:
            print("\nRestarting configuration...")
            return get_user_inputs()
        print("Please enter 'yes' or 'no'.")

    return folders, output_folder, target_masses, mass_tolerance, control_injection, alpha


def read_dat_file(filepath):
    """
    Read one UniDec DAT file into a DataFrame.
    """
    try:
        return pd.read_csv(filepath, sep=r"\s+", header=None, encoding="utf-8")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def get_peak_measurements(df, filepath):
    """
    Return numeric peakparam measurements needed for matching and statistics.
    """
    if df.shape[1] < REQUIRED_COLUMN_COUNT:
        print(
            f"Warning: Expected at least {REQUIRED_COLUMN_COUNT} columns in "
            f"{filepath}; found {df.shape[1]}"
        )
        return pd.DataFrame(columns=MEASUREMENT_COLUMNS.keys())

    measurements = pd.DataFrame({
        metric: pd.to_numeric(df.iloc[:, column_index], errors="coerce")
        for metric, column_index in MEASUREMENT_COLUMNS.items()
    })
    measurements = measurements.dropna(subset=["Mass Centroid"])

    if measurements.empty:
        print(f"Warning: No numeric mass centroid values found in {filepath}")

    return measurements


def extract_injection_number(filename):
    """
    Extract the injection number from names like '*merged_21_peakparam.dat'.
    """
    match = re.search(r"merged_(\d+)_", filename)
    return int(match.group(1)) if match else None


def match_centroid_to_target(mass_centroid, target_masses, tolerance):
    """
    Return the closest target mass to this centroid, or (None, None).
    """
    matched_target_mass = None
    matched_name = None
    min_diff = float("inf")

    for name, target_mass in target_masses.items():
        diff = abs(mass_centroid - target_mass)
        if diff < min_diff and diff <= tolerance:
            min_diff = diff
            matched_target_mass = target_mass
            matched_name = name

    return matched_target_mass, matched_name


def collect_matched_measurements(folders, target_masses, mass_tolerance):
    """
    Collect one closest matched peak per target mass from each DAT file.

    Area and height ratios are normalized to the total area/height of the
    matched target peaks in the same DAT file.
    """
    records = []

    for folder_index, folder in enumerate(folders, start=1):
        folder_path = Path(folder)

        if not folder_path.exists():
            print(f"Warning: Folder does not exist: {folder}")
            continue

        dat_files = sorted(folder_path.glob("*_peakparam.dat"))
        print(f"\nProcessing folder {folder_index}: Found {len(dat_files)} DAT files")

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

            closest_matches = {}
            for _, peak in measurements.iterrows():
                mass_centroid = float(peak["Mass Centroid"])
                target_mass, mass_name = match_centroid_to_target(
                    mass_centroid, target_masses, mass_tolerance
                )

                if target_mass is None:
                    continue

                candidate_diff = abs(mass_centroid - target_mass)
                current = closest_matches.get(mass_name)
                if current is None or candidate_diff < current["match_error_da"]:
                    closest_matches[mass_name] = {
                        "folder_index": folder_index,
                        "folder": str(folder_path),
                        "file": str(dat_file),
                        "file_name": dat_file.name,
                        "injection": injection,
                        "peak_name": mass_name,
                        "target_mass_da": target_mass,
                        "match_error_da": candidate_diff,
                        "Mass Centroid": mass_centroid,
                        "Area": float(peak["Area"]),
                        "Height": float(peak["Height"]),
                        "FWHM": float(peak["FWHM"]),
                    }

            if not closest_matches:
                continue

            total_area = sum(
                match["Area"] for match in closest_matches.values()
                if pd.notna(match["Area"])
            )
            total_height = sum(
                match["Height"] for match in closest_matches.values()
                if pd.notna(match["Height"])
            )

            if total_area <= 0:
                print(f"Warning: Matched peak areas sum to zero or less in {dat_file}")
            if total_height <= 0:
                print(f"Warning: Matched peak heights sum to zero or less in {dat_file}")

            for match in closest_matches.values():
                match["Area Ratio"] = (
                    match["Area"] / total_area
                    if total_area > 0 and pd.notna(match["Area"])
                    else np.nan
                )
                match["Height Ratio"] = (
                    match["Height"] / total_height
                    if total_height > 0 and pd.notna(match["Height"])
                    else np.nan
                )
                records.append(match)

    return pd.DataFrame.from_records(records)


def get_metric_groups(metric_df):
    """
    Return sorted injection groups as {injection: numpy_array}.
    """
    groups = {}
    for injection, group in metric_df.groupby("injection"):
        values = group["value"].dropna().to_numpy(dtype=float)
        if values.size:
            groups[int(injection)] = values
    return dict(sorted(groups.items()))


def run_one_way_anova(groups):
    """
    Run one-way ANOVA across injection groups when enough data is available.
    """
    total_n = sum(len(values) for values in groups.values())
    residual_df = total_n - len(groups)

    if len(groups) < 2:
        return np.nan, np.nan, "Need at least two injection groups."
    if residual_df <= 0:
        return np.nan, np.nan, "Need replicate observations for residual variance."

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        result = stats.f_oneway(*groups.values())

    return float(result.statistic), float(result.pvalue), ""


def run_dunnett(groups, control_injection, alpha):
    """
    Compare every injection group to the selected control injection using Dunnett's test.
    """
    if control_injection not in groups:
        return [], f"Control injection {control_injection} was not found."

    comparison_injections = [
        injection for injection in groups.keys()
        if injection != control_injection and len(groups[injection]) > 0
    ]
    if not comparison_injections:
        return [], "No non-control injection groups were available."

    total_n = sum(len(groups[injection]) for injection in comparison_injections)
    total_n += len(groups[control_injection])
    residual_df = total_n - (len(comparison_injections) + 1)
    if residual_df <= 0:
        return [], "Need replicate observations for Dunnett residual variance."

    try:
        result = stats.dunnett(
            *(groups[injection] for injection in comparison_injections),
            control=groups[control_injection],
            alternative="two-sided",
        )
        ci = result.confidence_interval(confidence_level=1 - alpha)
    except Exception as e:
        return [], f"Dunnett test failed: {e}"

    rows = []
    control_mean = float(np.mean(groups[control_injection]))

    for idx, injection in enumerate(comparison_injections):
        comparison_mean = float(np.mean(groups[injection]))
        p_value = float(result.pvalue[idx])
        rows.append({
            "control_injection": control_injection,
            "comparison_injection": injection,
            "control_n": int(len(groups[control_injection])),
            "comparison_n": int(len(groups[injection])),
            "control_mean": control_mean,
            "comparison_mean": comparison_mean,
            "mean_difference": comparison_mean - control_mean,
            "dunnett_statistic": float(result.statistic[idx]),
            "dunnett_p": p_value,
            "ci_low": float(ci.low[idx]),
            "ci_high": float(ci.high[idx]),
            "alpha": alpha,
            "significant": p_value < alpha,
            "note": "",
        })

    return rows, ""


def analyze_measurements(raw_measurements, target_masses, control_injection, alpha):
    """
    Run ANOVA and Dunnett tests for each target peak and measurement metric.
    """
    anova_rows = []
    dunnett_rows = []

    for peak_name in target_masses.keys():
        peak_data = raw_measurements[raw_measurements["peak_name"] == peak_name]

        for metric in STAT_METRICS:
            metric_df = peak_data[["injection", metric]].rename(columns={metric: "value"})
            metric_df = metric_df.dropna(subset=["value"])
            groups = get_metric_groups(metric_df)

            anova_f, anova_p, anova_note = run_one_way_anova(groups)
            anova_rows.append({
                "peak_name": peak_name,
                "metric": metric,
                "groups_tested": len(groups),
                "total_n": int(sum(len(values) for values in groups.values())),
                "anova_f": anova_f,
                "anova_p": anova_p,
                "note": anova_note,
            })

            comparison_rows, dunnett_note = run_dunnett(groups, control_injection, alpha)
            if comparison_rows:
                for row in comparison_rows:
                    row.update({
                        "peak_name": peak_name,
                        "metric": metric,
                        "anova_f": anova_f,
                        "anova_p": anova_p,
                    })
                    dunnett_rows.append(row)
            elif dunnett_note:
                dunnett_rows.append({
                    "peak_name": peak_name,
                    "metric": metric,
                    "control_injection": control_injection,
                    "comparison_injection": np.nan,
                    "control_n": int(len(groups.get(control_injection, []))),
                    "comparison_n": 0,
                    "control_mean": np.nan,
                    "comparison_mean": np.nan,
                    "mean_difference": np.nan,
                    "anova_f": anova_f,
                    "anova_p": anova_p,
                    "dunnett_statistic": np.nan,
                    "dunnett_p": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "alpha": alpha,
                    "significant": False,
                    "note": dunnett_note,
                })

    return pd.DataFrame(anova_rows), pd.DataFrame(dunnett_rows)


def save_results(raw_measurements, anova_results, dunnett_results, output_folder):
    """
    Save raw measurements and statistical results as CSV files.
    """
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    raw_file = output_path / "dunnett_anova_raw_measurements.csv"
    anova_file = output_path / "dunnett_anova_summary.csv"
    dunnett_file = output_path / "dunnett_anova_posthoc.csv"

    raw_measurements.to_csv(raw_file, index=False)
    anova_results.to_csv(anova_file, index=False)
    dunnett_results.to_csv(dunnett_file, index=False)

    return raw_file, anova_file, dunnett_file


def print_summary(raw_measurements, anova_results, dunnett_results, output_files):
    """
    Print a concise analysis summary.
    """
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"Matched peak measurements: {len(raw_measurements)}")

    significant = dunnett_results[dunnett_results["significant"] == True]
    print(f"Significant Dunnett comparisons: {len(significant)}")

    if not significant.empty:
        print("\nSignificant comparisons:")
        display_columns = [
            "peak_name", "metric", "comparison_injection",
            "control_injection", "mean_difference", "dunnett_p",
        ]
        print(significant[display_columns].to_string(index=False))

    print("\nANOVA rows with notes:")
    noted_anova = anova_results[anova_results["note"].fillna("") != ""]
    if noted_anova.empty:
        print("  None")
    else:
        print(noted_anova[["peak_name", "metric", "note"]].to_string(index=False))

    print("\nOutput files:")
    for output_file in output_files:
        print(f"  {output_file}")


def main():
    """
    Run the interactive statistics workflow.
    """
    folders, output_folder, target_masses, mass_tolerance, control_injection, alpha = get_user_inputs()

    print("\n" + "=" * 80)
    print("Starting DAT file statistical analysis...")
    print("=" * 80)

    raw_measurements = collect_matched_measurements(folders, target_masses, mass_tolerance)
    if raw_measurements.empty:
        print("\nNo matched peak measurements were found. Please check folder paths, file names, and target masses.")
        return

    anova_results, dunnett_results = analyze_measurements(
        raw_measurements,
        target_masses,
        control_injection,
        alpha,
    )
    output_files = save_results(
        raw_measurements,
        anova_results,
        dunnett_results,
        output_folder,
    )
    print_summary(raw_measurements, anova_results, dunnett_results, output_files)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
