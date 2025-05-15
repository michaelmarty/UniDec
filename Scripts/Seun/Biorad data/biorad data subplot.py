import matplotlib
matplotlib.use("TkAgg")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to create folder if it doesn't exist
def ensure_folder_exists(path):
    folder = os.path.dirname(path)
    if not os.path.exists(folder):
        os.makedirs(folder)

# Function to read CSV files and plot chromatograms safely
def plot_chromatograms(file_paths, legend_names, x_label=None, y_label=None):
    num_files = len(file_paths)
    num_subplots = 3  # Set how many plots per figure

    colors = ['black', 'gray', 'green', 'red', 'brown', 'limegreen', 'indianred', 'blue', 'olive', 'teal']

    for i in range(0, num_files, num_subplots):
        plt.figure(figsize=(10, 8))

        for j in range(num_subplots):
            idx = i + j
            if idx >= len(file_paths):
                break  # No more files to plot

            try:
                # Read the CSV
                data = pd.read_csv(file_paths[idx], header=1)

                volume_column_name = data.columns[0]
                uv_absorbance_column_name = data.columns[1]

                volume = data[volume_column_name]
                uv_absorbance = data[uv_absorbance_column_name]
                uv_absorbance /= uv_absorbance.max()  # Normalize

                # Safe label assignment
                label = legend_names[idx] if idx < len(legend_names) else f"File {idx+1}"
                color = colors[idx % len(colors)]

                plt.subplot(num_subplots, 1, j + 1)
                plt.plot(volume, uv_absorbance, color=color, label=label, linewidth=2)
                plt.xlabel(x_label if x_label else volume_column_name)
                plt.ylabel(y_label if y_label else uv_absorbance_column_name)
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.ylim(0, 1)
                plt.xlim(6, 24)
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.title(label)

            except Exception as e:
                print(f"Error plotting file {file_paths[idx]}: {e}")

        plt.tight_layout()

        # Safe saving
        save_path = f"Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\LC\\Chromatogram\\20250430_Nanodisc_chromatograms_subplot_{i//num_subplots + 1}.png"
        ensure_folder_exists(save_path)
        try:
            plt.savefig(save_path, dpi=600)
            print(f"Saved: {save_path}")
        except Exception as e:
            print(f"Error saving figure {save_path}: {e}")

        plt.close()

# === Define your files and legends ===

file_paths = [
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250419_OF_MSP1D1 POPC 75_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\control data\\20250426_OF_(+)MSP1E3D1 POPC 140_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\D1_sf GFP Data\\20250313_OF_(+)MSP1D1 sf GFP POPC 65_1 Rep2.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\D1_mCherry Data\\20250213_OF_(+)MSP1D1 mCherry POPC 65_1 Rep2.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\20250203\\20250203_OF_(+)MSP1D1 Halo POPC 80_1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\E3_sf GFP data\\20250426_OF_(+)MSP1E3D1_sfGFP POPC 160_1 rep 1.csv",
    "Z:\\Group Share\\Seun\\Current data\\MSP Engineering\\Biorad data\\E3_mCherry Data\\20250214_OF_(+)MSP1E3D1 mCherry POPC 130_1 Rep2.csv",
]

legend_names = [
    '(+)MSP1D1 POPC',
    '(+)MSP1E3D1 POPC',
    '(+)MSP1D1 sf GFP POPC',
    '(+)MSP1D1 mCherry POPC',
    '(+)MSP1D1 Halo POPC',
    '(+)MSP1E3D1 sf GFP POPC',
    '(+)MSP1E3D1 mCherry POPC'
]

custom_x_label = 'Retention Volume (ml)'
custom_y_label = 'A280 nm (mAU)'

# === Run the function ===

plot_chromatograms(file_paths, legend_names, x_label=custom_x_label, y_label=custom_y_label)
