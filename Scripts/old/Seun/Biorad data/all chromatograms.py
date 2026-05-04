import matplotlib
matplotlib.use("TkAgg")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to read CSV files and plot chromatograms with custom legend names and axis labels
def plot_chromatograms(file_paths, legend_names, x_label=None, y_label=None):

    # Define a color palette with enough variety
    colors = ['black', 'grey', 'green', 'red', 'brown', 'indianred', 'cyan', 'gray', 'olive', 'teal']

    # Create a single figure for all chromatograms
    plt.figure(figsize=(10, 6))

    for idx, file_path in enumerate(file_paths):
        # Read the CSV file
        data = pd.read_csv(file_path, header=1)

        # Use the actual column names from the CSV
        volume_column_name = data.columns[0]  # The first column name
        uv_absorbance_column_name = data.columns[1]  # The second column name

        # Extract the data
        volume = data[volume_column_name]
        uv_absorbance = data[uv_absorbance_column_name]

        # Assign a label and color based on the file index
        label = legend_names[idx]
        color = colors[idx % len(colors)]  # Ensure color rotation if more files than colors

        # Plot the data on the same axes
        plt.plot(volume, uv_absorbance, color=color, label=label, linewidth=2)

    # Customize the plot
    plt.xlabel(x_label if x_label else volume_column_name)
    plt.ylabel(y_label if y_label else uv_absorbance_column_name)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylim(0, 400)
    plt.xlim(6, 24)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Adjust layout for tightness
    plt.tight_layout()

    # Save the plot
    save_path = f"Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\LC\\Chromatogram\\202504224_Nanodisc chromatograms_overlaynorm.png"
    plt.savefig(save_path, dpi=1000)
    plt.show()
    plt.close()  # Close the figure to avoid displaying it again

# Replace the list with the actual paths to your CSV files
file_paths = [
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250419_OF_MSP1D1 POPC 75_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\EC (+)MSP1E3D1 POPC 130_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\D1_sf GFP Data\\20250313_OF_(+)MSP1D1 sf GFP POPC 65_1 Rep2.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\D1_mCherry Data\\20250213_OF_(+)MSP1D1 mCherry POPC 65_1 Rep2.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\20250203\\20250203_OF_(+)MSP1D1 Halo POPC 80_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\E3_mCherry Data\\20250214_OF_(+)MSP1E3D1 mCherry POPC 130_1 Rep2.csv",
]

# Provide custom legend names corresponding to each file path
legend_names = [
    '(+)MSP1D1 POPC',
    '(+)MSP1E3D1 POPC',
    '(+)MSP1D1 sf GFP POPC',
    '(+)MSP1D1 mCherry POPC',
    '(+)MSP1D1 Halo POPC',
    '(+)MSP1E3D1 mCherry POPC'
]

# Specify custom axis labels if desired
custom_x_label = 'Retention Volume (ml)'
custom_y_label = 'A280 nm (mAU)'

plot_chromatograms(file_paths, legend_names, x_label=custom_x_label, y_label=custom_y_label)
