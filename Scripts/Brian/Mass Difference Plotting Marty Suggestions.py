import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path


PEAK_REFERENCE_FILE_COUNT = 3
TARGET_INJECTION_REFERENCE = 21


def clean_path(path_string):
    """
    Clean a path string by removing quotes and extra whitespace
    Handles both single and double quotes from Windows copy-paste
    """
    # Strip whitespace
    path_string = path_string.strip()

    # Remove surrounding quotes (both single and double)
    if (path_string.startswith('"') and path_string.endswith('"')) or \
            (path_string.startswith("'") and path_string.endswith("'")):
        path_string = path_string[1:-1]

    # Strip again in case there was whitespace inside quotes
    path_string = path_string.strip()

    return path_string


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
    Get all user inputs for the analysis
    """
    print("=" * 80)
    print("DAT FILE MASS ANALYSIS")
    print("=" * 80)

    # Get number of replicates
    while True:
        try:
            num_reps = int(input("\nEnter the number of replicate folders: "))
            if num_reps > 0:
                break
            else:
                print("Please enter a positive number.")
        except ValueError:
            print("Please enter a valid integer.")

    # Get folder paths
    folders = prompt_for_paths(num_reps, "replicate folders")

    # Get peak text file paths used for the horizontal mean/SEM reference line
    peak_reference_files = prompt_for_paths(
        PEAK_REFERENCE_FILE_COUNT, "reference peak files"
    )

    # Get output folder
    print("\nEnter the output folder path:")
    print("(Tip: You can copy-paste paths from Windows Explorer - quotes will be handled automatically)")
    output_folder = input("  Output folder: ")
    output_folder = clean_path(output_folder)
    if not output_folder:
        output_folder = "."  # Current directory as default
        print(f"  Using current directory: {output_folder}")

    # Get number of theoretical masses
    while True:
        try:
            num_masses = int(input("\nEnter the number of theoretical masses to track: "))
            if num_masses > 0:
                break
            else:
                print("Please enter a positive number.")
        except ValueError:
            print("Please enter a valid integer.")

    # Get theoretical masses with names
    theoretical_masses = {}
    print(f"\nEnter {num_masses} theoretical mass values with their names:")
    for i in range(num_masses):
        print(f"\n  Mass {i + 1}:")
        name = input("    Name (e.g., 'Cytochrome C'): ").strip()
        while True:
            try:
                mass_value = float(input("    Theoretical mass (Da): "))
                if mass_value > 0:
                    theoretical_masses[name] = mass_value
                    break
                else:
                    print("    Please enter a positive mass value.")
            except ValueError:
                print("    Please enter a valid number.")

    # Get mass tolerance
    while True:
        try:
            mass_tolerance = float(input("\nEnter mass tolerance for matching (Da) [default: 50]: ").strip() or "50")
            if mass_tolerance > 0:
                break
            else:
                print("Please enter a positive tolerance value.")
        except ValueError:
            print("Please enter a valid number.")

    # Summary
    print("\n" + "=" * 80)
    print("CONFIGURATION SUMMARY")
    print("=" * 80)
    print(f"\nNumber of replicates: {num_reps}")
    print("\nReplicate folders:")
    for i, folder in enumerate(folders, 1):
        print(f"  {i}. {folder}")
    print("\nReference peak files:")
    for i, peak_file in enumerate(peak_reference_files, 1):
        print(f"  {i}. {peak_file}")
    print(f"\nOutput folder: {output_folder}")
    print(f"\nMass tolerance: {mass_tolerance} Da")
    print("\nTheoretical masses:")
    for name, mass in theoretical_masses.items():
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

    return folders, peak_reference_files, output_folder, theoretical_masses, mass_tolerance


def read_dat_file(filepath):
    """
    Read a DAT file and return a pandas DataFrame
    """
    try:
        df = pd.read_csv(filepath, sep=r'\s+', encoding='utf-8')
        return df
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def read_peak_masses(filepath):
    """
    Read a UTF-8 comma-separated UniDec peaks text file and return first-column masses.
    """
    masses = []
    try:
        with open(filepath, "r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            for row_num, row in enumerate(reader, start=1):
                if not row or not row[0].strip():
                    continue
                try:
                    masses.append(float(row[0].strip()))
                except ValueError:
                    print(
                        f"Warning: Could not parse first-column mass in "
                        f"{filepath} line {row_num}: {row[0]}"
                    )
    except Exception as e:
        print(f"Error reading reference peak file {filepath}: {e}")

    return masses


def match_mass_to_theoretical(measured_mass, theoretical_masses, tolerance=50):
    """
    Match a measured mass to the closest theoretical mass within tolerance
    Returns the theoretical mass value and the name, or (None, None) if no match
    """
    min_diff = float('inf')
    matched_theoretical = None
    matched_name = None

    for name, theo_mass in theoretical_masses.items():
        diff = abs(measured_mass - theo_mass)
        if diff < min_diff and diff <= tolerance:
            min_diff = diff
            matched_theoretical = theo_mass
            matched_name = name

    return matched_theoretical, matched_name


def match_peak_masses_to_theoretical(peak_masses, theoretical_masses, tolerance=50):
    """
    Match each theoretical protein to the closest peak mass within tolerance.
    """
    matches = {}
    for name, theo_mass in theoretical_masses.items():
        best_mass = None
        best_diff = float("inf")
        for peak_mass in peak_masses:
            diff = abs(peak_mass - theo_mass)
            if diff < best_diff and diff <= tolerance:
                best_diff = diff
                best_mass = peak_mass

        if best_mass is not None:
            matches[name] = best_mass

    return matches


def extract_injection_number(filename):
    """
    Extract the injection number from the filename
    Example: ...merged_21_peakparam.dat -> 21
    """
    try:
        # Split by 'merged_' and take the part after it
        parts = filename.split('merged_')
        if len(parts) > 1:
            # Get the number before '_peakparam.dat'
            num_str = parts[1].split('_')[0]
            return int(num_str)
    except:
        pass
    return None


def process_all_data(folders, theoretical_masses, mass_tolerance=50):
    """
    Process all DAT files from all folders
    Returns a dictionary with structure: {injection_num: {rep_num: {mass_name: measured_mass}}}
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

            if 'Mass' not in df.columns:
                print(f"Warning: No 'Mass' column found in {dat_file}")
                continue

            # Initialize data structure
            if inj_num not in all_data:
                all_data[inj_num] = {}

            if rep_idx not in all_data[inj_num]:
                all_data[inj_num][rep_idx] = {}

            # Match each measured mass to theoretical masses. If multiple peaks match the
            # same protein in a file, keep the closest measured mass.
            measured_masses = pd.to_numeric(df['Mass'], errors='coerce').dropna()
            for measured_mass in measured_masses:
                theo_mass, mass_name = match_mass_to_theoretical(
                    measured_mass, theoretical_masses, mass_tolerance
                )

                if theo_mass is not None:
                    current_mass = all_data[inj_num][rep_idx].get(mass_name)
                    if current_mass is None or abs(measured_mass - theo_mass) < abs(current_mass - theo_mass):
                        all_data[inj_num][rep_idx][mass_name] = float(measured_mass)

    return all_data


def calculate_reference_peak_statistics(peak_reference_files, theoretical_masses, mass_tolerance=50):
    """
    Calculate mean mass and SEM for each protein from the reference peak files.
    """
    reference_values = {mass_name: [] for mass_name in theoretical_masses.keys()}

    for peak_file in peak_reference_files:
        peak_path = Path(peak_file)
        if not peak_path.exists():
            print(f"Warning: Reference peak file does not exist: {peak_file}")
            continue

        peak_masses = read_peak_masses(peak_path)
        if not peak_masses:
            print(f"Warning: No parseable masses found in reference peak file: {peak_file}")
            continue

        matches = match_peak_masses_to_theoretical(
            peak_masses, theoretical_masses, mass_tolerance
        )
        for mass_name, measured_mass in matches.items():
            reference_values[mass_name].append(measured_mass)

        missing_names = sorted(set(theoretical_masses.keys()) - set(matches.keys()))
        for mass_name in missing_names:
            print(
                f"Warning: No reference peak matched {mass_name} within "
                f"{mass_tolerance} Da in {peak_file}"
            )

    reference_stats = {}
    for mass_name, values in reference_values.items():
        if values:
            mean_val = float(np.mean(values))
            se_val = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
            reference_stats[mass_name] = {
                'mean': mean_val,
                'se': se_val,
                'n': len(values),
                'values': values,
            }
        else:
            reference_stats[mass_name] = {
                'mean': None,
                'se': None,
                'n': 0,
                'values': [],
            }

    return reference_stats


def calculate_statistics(all_data, theoretical_masses, num_replicates):
    """
    Calculate mean and standard error for each injection number and mass
    Returns DataFrames for plotting
    """
    results = {mass_name: {'injections': [], 'mean': [], 'se': []}
               for mass_name in theoretical_masses.keys()}

    # Get sorted injection numbers
    injection_numbers = sorted(all_data.keys())

    for inj_num in injection_numbers:
        for mass_name in theoretical_masses.keys():
            # Collect values from all replicates
            values = []
            for rep_num in range(1, num_replicates + 1):
                if rep_num in all_data[inj_num]:
                    if mass_name in all_data[inj_num][rep_num]:
                        values.append(all_data[inj_num][rep_num][mass_name])

            if values:
                mean_val = np.mean(values)
                # Standard error = std / sqrt(n)
                se_val = np.std(values, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0

                results[mass_name]['injections'].append(inj_num)
                results[mass_name]['mean'].append(mean_val)
                results[mass_name]['se'].append(se_val)

    return results


def set_mass_axis_limits(ax, data, reference_stats):
    """
    Set y-axis limits around measured masses and reference mean +/- SEM.
    """
    y_values = []

    if data['mean']:
        means = np.asarray(data['mean'], dtype=float)
        ses = np.asarray(data['se'], dtype=float)
        y_values.extend((means - ses).tolist())
        y_values.extend((means + ses).tolist())

    if reference_stats and reference_stats['mean'] is not None:
        ref_mean = reference_stats['mean']
        ref_se = reference_stats['se'] or 0.0
        y_values.extend([ref_mean - ref_se, ref_mean + ref_se])

    if not y_values:
        return

    ymin = min(y_values)
    ymax = max(y_values)
    y_span = ymax - ymin
    padding = max(y_span * 0.08, max(abs(ymax), 1) * 0.00005, 1)
    ax.set_ylim(ymin - padding, ymax + padding)


def get_injection_reference_stats(data, injection_number):
    """
    Get plotted mean and SEM values for one injection number.
    """
    for injection, mean_val, se_val in zip(data['injections'], data['mean'], data['se']):
        if injection == injection_number:
            return float(mean_val), float(se_val)
    return None, None


def add_reference_mean_sem_line(ax, ref_mean, ref_se, color, label):
    """
    Draw one dashed reference line centered on the mean with +/- SEM thickness.
    """
    if ref_se <= 0:
        ax.axhline(
            ref_mean,
            color=color,
            linestyle='--',
            linewidth=2.5,
            label=label,
            zorder=1
        )
        return

    dash_width = 0.045
    gap_width = dash_width
    x_start = 0.0
    first_dash = True
    while x_start < 1.0:
        width = min(dash_width, 1.0 - x_start)
        ax.add_patch(
            Rectangle(
                (x_start, ref_mean - ref_se),
                width,
                2 * ref_se,
                transform=ax.get_yaxis_transform(),
                facecolor=color,
                edgecolor='none',
                label=label if first_dash else None,
                zorder=1,
                clip_on=True
            )
        )
        first_dash = False
        x_start += dash_width + gap_width


def plot_results(results, theoretical_masses, reference_stats, output_folder):
    """
    Create stacked line plots with error bars
    """
    # Set font to Arial
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 18
    # Enable editable text in PDF
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    num_plots = len(theoretical_masses)

    # Calculate figure size for 16:9 aspect ratio
    # Width = 16 inches, Height = 9 inches
    fig, axes = plt.subplots(num_plots, 1, figsize=(16, 9))

    # Handle case where there's only one plot
    if num_plots == 1:
        axes = [axes]

    # Add main title with larger font
    fig.suptitle('Measured Mass vs Number of Injections',
                 fontsize=24, fontweight='bold', y=0.98)

    # Use a single consistent color for the marker points and lines across all subplots
    point_color = plt.get_cmap('tab10')(0)  # first color from the tab10 colormap
    reference_color = plt.get_cmap('tab10')(3)
    injection_reference_color = plt.get_cmap('tab10')(2)

    for idx, (mass_name, theo_mass) in enumerate(theoretical_masses.items()):
        ax = axes[idx]
        data = results[mass_name]
        reference = reference_stats.get(mass_name)
        show_legend = False

        if reference and reference['mean'] is not None:
            ref_mean = reference['mean']
            ref_se = reference['se'] or 0.0
            add_reference_mean_sem_line(
                ax,
                ref_mean,
                ref_se,
                reference_color,
                f'Reference mean +/- SEM (n={reference["n"]})'
            )
            show_legend = True

        inj_mean, inj_se = get_injection_reference_stats(data, TARGET_INJECTION_REFERENCE)
        if inj_mean is not None and inj_se > 0:
            ax.axhspan(
                inj_mean - inj_se, inj_mean + inj_se,
                facecolor=injection_reference_color, alpha=0.45, linewidth=0,
                label=f'{TARGET_INJECTION_REFERENCE} injections mean +/- SEM',
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

        set_mass_axis_limits(ax, data, reference)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Only add x-axis label to bottom subplot with larger font
        if idx == num_plots - 1:
            ax.set_xlabel('Number of Injections', fontsize=24, fontweight='bold')

        # Set subtitle with protein name
        ax.set_title(f'{mass_name} (Theoretical Mass: {theo_mass:.2f} Da)',
                     fontsize=18, fontweight='bold', pad=10)

        if show_legend:
            ax.legend(frameon=False, fontsize=12, loc='best')

    # Add unified y-axis label with larger font
    fig.text(0.04, 0.5, 'Mass (Da)',
             va='center', rotation='vertical', fontsize=24, fontweight='bold')

    plt.tight_layout(rect=(0.06, 0, 1, 0.97))  # Adjust layout to make room for y-label and title

    # Create output folder if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save as PDF with editable text
    output_file = output_path / 'mass_analysis.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight',
                transparent=True, format='pdf')
    plt.show()

    print(f"\nPlot saved as '{output_file}'")
    print("Format: PDF (300 DPI, transparent, 16:9 aspect ratio, editable text)")


def create_summary_table(results, theoretical_masses, reference_stats):
    """
    Create a summary table of the results
    """
    print("\n" + "=" * 80)
    print("SUMMARY TABLE: Mean Measured Masses (Da)")
    print("=" * 80)

    for mass_name, theo_mass in theoretical_masses.items():
        print(f"\n{mass_name} (Theoretical: {theo_mass:.2f} Da)")
        print("-" * 60)
        print(f"{'Injections':<12} {'Mean Mass (Da)':<18} {'Std Error (Da)':<18}")
        print("-" * 60)

        data = results[mass_name]
        for inj, mean_val, se_val in zip(data['injections'], data['mean'], data['se']):
            print(f"{inj:<12} {mean_val:<18.4f} {se_val:<18.4f}")

        reference = reference_stats.get(mass_name)
        if reference and reference['mean'] is not None:
            print(
                f"Reference peak files: mean={reference['mean']:.4f} Da, "
                f"SEM={reference['se']:.4f} Da, n={reference['n']}"
            )
        else:
            print("Reference peak files: no matched masses")


# Main execution
if __name__ == "__main__":
    # Get user inputs
    folders, peak_reference_files, output_folder, theoretical_masses, mass_tolerance = get_user_inputs()

    print("\n" + "=" * 80)
    print("Starting DAT file analysis...")
    print("=" * 80)

    # Process all data
    all_data = process_all_data(folders, theoretical_masses, mass_tolerance)

    if not all_data:
        print("\nNo data was processed. Please check folder paths and file names.")
    else:
        # Calculate statistics
        results = calculate_statistics(all_data, theoretical_masses, len(folders))
        reference_stats = calculate_reference_peak_statistics(
            peak_reference_files, theoretical_masses, mass_tolerance
        )

        # Create summary table
        create_summary_table(results, theoretical_masses, reference_stats)

        # Plot results
        plot_results(results, theoretical_masses, reference_stats, output_folder)

        print("\n" + "=" * 80)
        print("Analysis complete!")
        print("=" * 80)
