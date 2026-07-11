"""
Read UniDec *_peakparam.dat files from replicate folders, match mass centroids
to user-provided target masses, and plot mean height fractions +/- SEM by injection.
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from pathlib import Path


DEFAULT_MASS_TOLERANCE = 50.0
HIGHLIGHT_INJECTION = 21
MASS_CENTROID_COLUMN_INDEX = 6
MASS_CENTROID_COLUMN_NUMBER = MASS_CENTROID_COLUMN_INDEX + 1
HEIGHT_COLUMN_INDEX = 4
HEIGHT_COLUMN_NUMBER = HEIGHT_COLUMN_INDEX + 1
Y_TICK_INTERVALS_PER_SIDE = 2


def clean_path(path_string):
    """
    Normalize a path pasted from Windows Explorer or the terminal.
    """
    path_string = path_string.strip()

    if (path_string.startswith('"') and path_string.endswith('"')) or \
            (path_string.startswith("'") and path_string.endswith("'")):
        path_string = path_string[1:-1]

    return path_string.strip()


def prompt_for_positive_int(prompt):
    """
    Prompt until the user enters a positive integer.
    """
    while True:
        try:
            value = int(input(prompt))
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
    Prompt for a fixed number of paths.
    """
    paths = []
    print(f"\nEnter the paths for {count} {prompt_name}:")
    print("(Tip: You can copy-paste paths from Windows Explorer - quotes will be handled automatically)")
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
    Get all user inputs needed for the analysis.
    """
    print("=" * 80)
    print("DAT FILE HEIGHT RATIO ANALYSIS")
    print("=" * 80)

    num_reps = prompt_for_positive_int("\nEnter the number of replicate folders: ")
    folders = prompt_for_paths(num_reps, "replicate folders")

    print("\nEnter the output folder path:")
    print("(Tip: You can copy-paste paths from Windows Explorer - quotes will be handled automatically)")
    output_folder = input("  Output folder: ")
    output_folder = clean_path(output_folder)
    if not output_folder:
        output_folder = "."  # Current directory as default
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

    print("\n" + "=" * 80)
    print("CONFIGURATION SUMMARY")
    print("=" * 80)
    print(f"\nNumber of replicates: {num_reps}")
    print("\nReplicate folders:")
    for i, folder in enumerate(folders, 1):
        print(f"  {i}. {folder}")
    print(f"\nOutput folder: {output_folder}")
    print(f"\nMass tolerance: {mass_tolerance} Da")
    print("\nPeak targets:")
    for name, mass in target_masses.items():
        print(f"  {name}: {mass:.2f} Da")

    # Confirmation
    while True:
        confirm = input("\nProceed with this configuration? (yes/no): ").strip().lower()
        if confirm in ['yes', 'y']:
            break
        elif confirm in ['no', 'n']:
            print("\nRestarting configuration...")
            return get_user_inputs()
        else:
            print("Please enter 'yes' or 'no'.")

    return folders, output_folder, target_masses, mass_tolerance


def standard_error(values):
    """
    Calculate SEM; one value has no spread, so SEM is zero.
    """
    return float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0


def read_dat_file(filepath):
    """
    Read one UniDec DAT file into a DataFrame.
    """
    try:
        return pd.read_csv(filepath, sep=r'\s+', encoding='utf-8')
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def get_peak_measurements(df, filepath):
    """
    Return numeric mass centroid and height values from the peakparam DAT columns.
    """
    required_column_count = max(MASS_CENTROID_COLUMN_NUMBER, HEIGHT_COLUMN_NUMBER)
    if df.shape[1] < required_column_count:
        print(
            f"Warning: Expected at least {required_column_count} columns for "
            f"mass centroid/height in {filepath}; found {df.shape[1]}"
        )
        return pd.DataFrame(columns=['mass_centroid', 'height'])

    measurements = pd.DataFrame({
        'mass_centroid': pd.to_numeric(df.iloc[:, MASS_CENTROID_COLUMN_INDEX], errors='coerce'),
        'height': pd.to_numeric(df.iloc[:, HEIGHT_COLUMN_INDEX], errors='coerce'),
    }).dropna()
    if measurements.empty:
        print(f"Warning: No numeric mass centroid/height pairs found in {filepath}")

    return measurements


def match_centroid_to_target(mass_centroid, target_masses, tolerance=50):
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


def extract_injection_number(filename):
    """
    Extract the injection number from names like '*merged_21_peakparam.dat'.
    Example: ...merged_21_peakparam.dat -> 21
    """
    match = re.search(r"merged_(\d+)_", filename)
    return int(match.group(1)) if match else None


def process_all_data(folders, target_masses, mass_tolerance=50):
    """
    Return nested data as {injection_num: {rep_num: {mass_name: height_fraction}}}.
    """
    all_data = {}

    for rep_idx, folder in enumerate(folders, start=1):
        folder_path = Path(folder)

        if not folder_path.exists():
            print(f"Warning: Folder does not exist: {folder}")
            continue

        # Find all DAT files in the folder
        dat_files = sorted(folder_path.glob("*_peakparam.dat"))

        print(f"\nProcessing Rep {rep_idx}: Found {len(dat_files)} DAT files")

        for dat_file in dat_files:
            # Extract injection number
            inj_num = extract_injection_number(dat_file.name)

            if inj_num is None:
                print(f"Warning: Could not extract injection number from {dat_file.name}")
                continue

            # Read the DAT file
            df = read_dat_file(dat_file)

            if df is None or df.empty:
                continue

            peak_measurements = get_peak_measurements(df, dat_file)
            if peak_measurements.empty:
                continue

            matched_heights = {}
            matched_diffs = {}

            # If multiple peaks match the same protein, keep the height from the closest centroid.
            for peak in peak_measurements.itertuples(index=False):
                mass_centroid = float(peak.mass_centroid)
                target_mass, mass_name = match_centroid_to_target(
                    mass_centroid, target_masses, mass_tolerance
                )

                if target_mass is not None:
                    current_diff = matched_diffs.get(mass_name, float("inf"))
                    candidate_diff = abs(mass_centroid - target_mass)
                    if candidate_diff < current_diff:
                        matched_heights[mass_name] = float(peak.height)
                        matched_diffs[mass_name] = candidate_diff

            if not matched_heights:
                continue

            total_height = sum(matched_heights.values())
            if total_height <= 0:
                print(f"Warning: Matched peak heights sum to zero or less in {dat_file}")
                continue

            rep_data = all_data.setdefault(inj_num, {}).setdefault(rep_idx, {})
            for mass_name, height in matched_heights.items():
                rep_data[mass_name] = height / total_height

    return all_data


def calculate_statistics(all_data, target_masses, num_replicates):
    """
    Calculate mean and standard error for each injection number and height fraction.
    """
    results = {mass_name: {'injections': [], 'mean': [], 'se': []}
               for mass_name in target_masses.keys()}

    for inj_num in sorted(all_data.keys()):
        for mass_name in target_masses.keys():
            values = []
            for rep_num in range(1, num_replicates + 1):
                rep_values = all_data[inj_num].get(rep_num, {})
                if mass_name in rep_values:
                    values.append(rep_values[mass_name])

            if values:
                mean_val = float(np.mean(values))
                se_val = standard_error(values)

                results[mass_name]['injections'].append(inj_num)
                results[mass_name]['mean'].append(mean_val)
                results[mass_name]['se'].append(se_val)

    return results


def get_nice_tick_step(raw_step):
    """
    Round a raw tick interval up to a readable 1, 2, 2.5, 5, or 10 multiple.
    """
    if raw_step <= 0:
        return 1.0

    magnitude = 10 ** np.floor(np.log10(raw_step))
    normalized_step = raw_step / magnitude

    for candidate in (1, 2, 2.5, 5, 10):
        if normalized_step <= candidate:
            return float(candidate * magnitude)

    return float(10 * magnitude)


def set_ratio_axis_limits(ax, data):
    """
    Scale the y-axis around measured height fractions and end on tick marks.
    """
    if data['mean']:
        means = np.asarray(data['mean'], dtype=float)
        ses = np.asarray(data['se'], dtype=float)
        y_values = np.concatenate((means - ses, means + ses))
        lower_value = float(np.min(y_values))
        upper_value = float(np.max(y_values))
        span = upper_value - lower_value
        padding = max(span * 0.08, 0.01)
        lower_limit = lower_value - padding
        upper_limit = upper_value + padding

        if lower_limit == upper_limit:
            lower_limit -= 0.01
            upper_limit += 0.01

        tick_step = get_nice_tick_step((upper_limit - lower_limit) / (Y_TICK_INTERVALS_PER_SIDE * 2))
        tick_start = np.floor(lower_limit / tick_step) * tick_step
        tick_end = np.ceil(upper_limit / tick_step) * tick_step
        ticks = np.arange(tick_start, tick_end + tick_step * 0.5, tick_step)
    else:
        ticks = np.linspace(0, 1, 5)

    ax.set_ylim(ticks[0], ticks[-1])
    ax.set_yticks(ticks)

    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))
    ax.yaxis.get_offset_text().set_visible(False)


def get_highlight_band(data, injection_number):
    """
    Get the plotted mean and SEM values for one injection number.
    """
    for injection, mean_val, se_val in zip(data['injections'], data['mean'], data['se']):
        if injection == injection_number:
            return float(mean_val), float(se_val)
    return None, None


def plot_results(results, target_masses, output_folder):
    """
    Create stacked line plots with error bars.
    """
    # Set font to Arial
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 18
    # Enable editable text in PDF
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    num_plots = len(target_masses)

    # Calculate figure size for 16:9 aspect ratio
    # Width = 16 inches, Height = 9 inches
    fig, axes = plt.subplots(num_plots, 1, figsize=(16, 9))

    # Handle case where there's only one plot
    if num_plots == 1:
        axes = [axes]

    # Add main title with larger font
    fig.suptitle('Measured Height Fractions vs Number of Injections',
                 fontsize=24, fontweight='bold', y=0.98)

    # Use a single consistent color for the marker points and lines across all subplots
    point_color = plt.get_cmap('tab10')(0)  # first color from the tab10 colormap
    highlight_color = plt.get_cmap('tab10')(2)

    for idx, mass_name in enumerate(target_masses.keys()):
        ax = axes[idx]
        data = results[mass_name]
        show_legend = False

        inj_mean, inj_se = get_highlight_band(data, HIGHLIGHT_INJECTION)
        if inj_mean is not None and inj_se > 0:
            ax.axhspan(
                inj_mean - inj_se, inj_mean + inj_se,
                facecolor=highlight_color, alpha=0.45, linewidth=0,
                label=f'{HIGHLIGHT_INJECTION} injections mean +/- SEM',
                zorder=1
            )
            show_legend = True

        if data['injections']:
            # Keep error bars black while using a single consistent color for markers and lines
            ax.errorbar(
                data['injections'], data['mean'], yerr=data['se'],
                marker='o', linestyle='-', linewidth=2, markersize=8,
                capsize=5, capthick=2, color=point_color,
                markerfacecolor=point_color, markeredgecolor=point_color,
                ecolor='black', elinewidth=2, zorder=2
            )  # Black error bars

            # Set x-axis to show integer values
            ax.set_xticks(data['injections'])

        set_ratio_axis_limits(ax, data)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Only add x-axis label to bottom subplot with larger font
        if idx == num_plots - 1:
            ax.set_xlabel('Number of Injections', fontsize=24, fontweight='bold')

        # Set subtitle with the user-designated peak name only.
        ax.set_title(mass_name, fontsize=18, fontweight='bold', pad=10)

        if show_legend:
            ax.legend(frameon=False, fontsize=12, loc='best')

    # Add unified y-axis label with larger font
    fig.text(0.04, 0.5, 'Height Fraction',
             va='center', rotation='vertical', fontsize=24, fontweight='bold')

    plt.tight_layout(rect=(0.06, 0, 1, 0.97))  # Adjust layout to make room for y-label and title

    # Create output folder if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save as PDF with editable text
    output_file = output_path / 'height_ratio_analysis.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight',
                transparent=True, format='pdf')
    plt.show()

    print(f"\nPlot saved as '{output_file}'")
    print("Format: PDF (300 DPI, transparent, 16:9 aspect ratio, editable text)")


def create_summary_table(results, target_masses):
    """
    Print a summary table of the results.
    """
    print("\n" + "=" * 80)
    print("SUMMARY TABLE: Mean Height Fractions")
    print("=" * 80)

    for mass_name in target_masses.keys():
        print(f"\n{mass_name}")
        print("-" * 60)
        print(f"{'Injections':<12} {'Mean Fraction':<20} {'Std Error':<18}")
        print("-" * 60)

        data = results[mass_name]
        for inj, mean_val, se_val in zip(data['injections'], data['mean'], data['se']):
            print(f"{inj:<12} {mean_val:<20.4f} {se_val:<18.4f}")


def main():
    """
    Run the interactive analysis workflow.
    """
    folders, output_folder, target_masses, mass_tolerance = get_user_inputs()

    print("\n" + "=" * 80)
    print("Starting DAT file analysis...")
    print("=" * 80)

    # Process all data
    all_data = process_all_data(folders, target_masses, mass_tolerance)

    if not all_data:
        print("\nNo data was processed. Please check folder paths and file names.")
    else:
        # Calculate statistics
        results = calculate_statistics(all_data, target_masses, len(folders))

        # Create summary table
        create_summary_table(results, target_masses)

        # Plot results
        plot_results(results, target_masses, output_folder)

        print("\n" + "=" * 80)
        print("Analysis complete!")
        print("=" * 80)


if __name__ == "__main__":
    main()
