import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os


def prompt_yes_no(prompt):
    while True:
        response = input(prompt).strip().lower()
        if response in ("y", "yes"):
            return True
        if response in ("n", "no"):
            return False
        print("Please enter y or n.")


def prompt_zoom_range():
    while True:
        try:
            zoom_start = float(input("Enter the start mass for zoomed x-axis (in kDa): "))
            zoom_end = float(input("Enter the end mass for zoomed x-axis (in kDa): "))
        except ValueError:
            print("Please enter numeric values for the start and end masses.")
            continue

        if zoom_end <= zoom_start:
            print("End mass must be greater than start mass.")
            continue

        return zoom_start, zoom_end


def normalize_intensity(intensity, dataset_number):
    max_intensity = np.max(intensity)
    if max_intensity <= 0:
        raise ValueError(f"Cannot normalize dataset {dataset_number}; maximum intensity is {max_intensity}.")
    return intensity / max_intensity


# Prompt the user for the number of datasets
num_datasets = int(input("Enter the number of datasets: "))

# Initialize lists to store data and labels
file_paths = []
data_list = []
legend_labels = []

# Prompt the user for each file path, load the data, and get the legend label
for i in range(num_datasets):
    file_path = input(f"Enter the file path for dataset {i + 1} (with quotes): ").strip('"').strip("'")
    file_paths.append(file_path)
    data_list.append(np.loadtxt(file_path))
    legend_label = input(f"Enter the legend label for dataset {i + 1}: ")
    legend_labels.append(legend_label)

normalize_data = prompt_yes_no("Do you want to normalize the data? (y/n): ")
stagger_plots = prompt_yes_no("Do you want to stagger the plots? (y/n): ")
use_zoomed_x_range = prompt_yes_no("Do you want to use a zoomed mass x-axis range? (y/n): ")
if use_zoomed_x_range:
    zoom_start, zoom_end = prompt_zoom_range()
else:
    zoom_start = None
    zoom_end = None

# Set font properties and styling
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28
plt.rcParams['xtick.major.size'] = 12
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 12
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
plt.rcParams['axes.labelpad'] = 18
plt.rcParams['axes.titlepad'] = 18
plt.rcParams['axes.linewidth'] = plt.rcParams['xtick.major.width']
plt.rcParams['pdf.fonttype'] = 42  # Ensure fonts are saved as editable text in PDFs

# Create a new figure
plt.figure(figsize=(16, 9))

all_intensity = np.concatenate([data[:, 1] for data in data_list])
if normalize_data:
    stagger_offset_step = 0.1
else:
    stagger_offset_step = 0.1 * np.max(all_intensity)

plotted_intensity_list = []

# Loop through each dataset and plot
for i, data in enumerate(data_list):
    # Assuming the data has two columns: mass and intensity
    mass = data[:, 0]
    intensity = data[:, 1]

    # Convert mass from Da to kDa
    mass = mass / 1000

    if normalize_data:
        intensity = normalize_intensity(intensity, i + 1)

    if stagger_plots:
        intensity = intensity + (i * stagger_offset_step)

    plotted_intensity_list.append(intensity)

    # Plot the data with the user-provided legend label
    plt.plot(mass, intensity, label=legend_labels[i])

# Set labels
plt.xlabel(r'Mass (kDa)')
y_label = 'Normalized Intensity' if normalize_data else 'Intensity'
if stagger_plots:
    y_label += ' + Offset'
plt.ylabel(y_label)

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limits based on the plotted intensity values
if stagger_plots:
    plotted_intensity = np.concatenate(plotted_intensity_list)
    ax.set_ylim(0, np.max(plotted_intensity))
elif normalize_data:
    ax.set_ylim(0, 1.0)
else:
    ax.set_ylim(0, np.max(all_intensity))

# Set x-axis limits based on the user-defined zoom range or the combined datasets
if use_zoomed_x_range:
    ax.set_xlim(zoom_start, zoom_end)
else:
    all_mass = np.concatenate([data[:, 0] for data in data_list])
    ax.set_xlim(np.min(all_mass) / 1000, np.max(all_mass) / 1000)

# Set x-axis to display numbers without scientific notation
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Add a legend
plt.legend()

# Automatically adjust the layout to prevent label cropping
plt.tight_layout()

# Save the plot as a 300 DPI PNG with a transparent background in the same folder as the first dataset
output_name = 'mass_distribution_comparison'
if normalize_data:
    output_name += '_normalized'
if stagger_plots:
    output_name += '_staggered'
if use_zoomed_x_range:
    output_name += '_zoomed'

output_path_png = os.path.join(os.path.dirname(file_paths[0]), f'{output_name}.png')
plt.savefig(output_path_png, dpi=300, transparent=True)

# Save the plot as a PDF in the same folder as the first dataset
output_path_pdf = os.path.join(os.path.dirname(file_paths[0]), f'{output_name}.pdf')
plt.savefig(output_path_pdf, dpi=300)

# Print the output paths
print(f"Plot saved as PNG at: {output_path_png}")
print(f"Plot saved as PDF at: {output_path_pdf}")

# Show the plot
plt.show()
