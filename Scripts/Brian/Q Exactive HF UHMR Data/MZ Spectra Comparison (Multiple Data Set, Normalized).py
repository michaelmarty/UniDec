import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Prompt the user for the number of datasets
num_datasets = int(input("Enter the number of datasets: "))

# Initialize lists to store data
file_paths = []
data_list = []

# Prompt the user for each file path
for i in range(num_datasets):
    file_path = input(f"Enter the file path for dataset {i + 1} (with quotes): ").strip('"').strip("'")
    file_paths.append(file_path)
    data_list.append(np.loadtxt(file_path))

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

# Plot each dataset
for i, data in enumerate(data_list):
    # Assuming the data has two columns: m/z and intensity
    mz = data[:, 0]
    intensity = data[:, 1]

    # Normalize the intensity
    intensity = intensity / np.max(intensity)

    # Plot the data with a unique color for each dataset
    plt.plot(mz, intensity, label=f'Dataset {i + 1}')

# Set labels
plt.xlabel(r'$m/z$')
plt.ylabel('Normalized Intensity')

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limits to 1.0
ax.set_ylim(0, 1.0)

# Set x-axis limits based on the combined dataset
all_mz = np.concatenate([data[:, 0] for data in data_list])
ax.set_xlim(np.min(all_mz), np.max(all_mz))

# Set x-axis to display numbers without scientific notation
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Add a legend
plt.legend()

# Automatically adjust the layout to prevent label cropping
plt.tight_layout()

# Save the plot as a 300 DPI PNG with a transparent background in the same folder as the first dataset
output_path_png = os.path.join(os.path.dirname(file_paths[0]), 'mass_spectrum_overlay.png')
plt.savefig(output_path_png, dpi=300, transparent=True)

# Save the plot as a PDF in the same folder as the first dataset
output_path_pdf = os.path.join(os.path.dirname(file_paths[0]), 'mass_spectrum_overlay.pdf')
plt.savefig(output_path_pdf, dpi=300)

# Print the output paths
print(f"Overlay plot saved as PNG at: {output_path_png}")
print(f"Overlay plot saved as PDF at: {output_path_pdf}")

# Show the plot
plt.show()