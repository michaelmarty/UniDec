import matplotlib
matplotlib.use("TkAgg")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to read CSV files and plot chromatograms as overlayed plots
def plot_chromatograms(file_paths, legend_names, x_label=None, y_label=None):

    plt.figure(figsize=(10, 6))  # Create one figure for overlayed plots

    # Define different linestyles for each file
    linestyles = ['-', '--', ':']  # Solid, dashed, dotted

    for idx, file_path in enumerate(file_paths):
        # Read the CSV file, specifying the row containing column names (0-indexed)
        data = pd.read_csv(file_path, header=1)

        # Identify column names
        volume_column_name = data.columns[0]
        uv_absorbance_column_name = data.columns[1]

        volume = data[volume_column_name]
        uv_absorbance = data[uv_absorbance_column_name]

        # Assign label, color, and linestyle
        label = legend_names[idx]
        color = "black"
        linestyle = linestyles[idx % len(linestyles)]

        # Plot the chromatogram
        plt.plot(volume, uv_absorbance, color=color, linestyle=linestyle, label=label, linewidth=2)

    # Labels and appearance
    plt.xlabel(x_label if x_label else volume_column_name)
    plt.ylabel(y_label if y_label else uv_absorbance_column_name)
    plt.ylim(0, 400)
    plt.xlim(6, 24)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Legend outside the plot
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()

    # Save and show
    save_path = "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\LC\\Chromatogram\\20250428_(+)MSP1E3D1 POPC_overlay.png"
    plt.savefig(save_path, dpi=1000)
    plt.show()
    plt.close()

# Replace the list with the actual paths to your CSV files
file_paths = [
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250426_OF_(+)MSP1E3D1 POPC 130_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250426_OF_(+)MSP1E3D1 POPC 140_1.csv",
    #"Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250419_OF_(+)MSP1D1 POPC 85_1.csv",
]

# Provide custom legend names corresponding to each file path
legend_names = [
    '(+)MSP1ED1 POPC 130_1',
    '(+)MSP1E3D1 POPC 140_1',
    #'(+)MSP1D1 POPC 85_1',
]

# Custom axis labels
custom_x_label = 'Retention Volume (ml)'
custom_y_label = 'A280 nm (mAU)'

plot_chromatograms(file_paths, legend_names, x_label=custom_x_label, y_label=custom_y_label)
