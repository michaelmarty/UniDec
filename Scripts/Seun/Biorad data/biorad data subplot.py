import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to read CSV files and plot chromatograms with custom legend names and axis labels
def plot_chromatograms(file_paths, legend_names, x_label=None, y_label=None):

    num_files = len(file_paths)
    num_subplots = 3  # Set the number of subplots per figure

    for i in range(0, num_files, num_subplots):
        plt.figure(figsize=(10, 6))

        # Loop through the next three files
        for j in range(num_subplots):
            idx = i + j
            if idx >= num_files:
                break  # Stop if there are no more files

            # Read the CSV file, specifying the row containing column names (0-indexed)
            data = pd.read_csv(file_paths[idx], header=1)

            # Update these variable names based on the actual column names in your CSV
            volume_column_name = data.columns[0]  # The first column name
            uv_absorbance_column_name = data.columns[1]  # The second column name

            # Assuming the CSV file has the correct columns
            volume = data[volume_column_name]
            uv_absorbance = data[uv_absorbance_column_name]

            # Use the provided legend name for each plot
            label = legend_names[idx]
            color = "green" if idx == 0 else "indianred"  # First file in black, others in green

            plt.subplot(num_subplots, 1, j + 1)
            plt.plot(volume, uv_absorbance, color=color, label=label, linewidth=2)

            # Add titles and labels
            plt.xlabel(x_label if x_label else volume_column_name)
            plt.ylabel(y_label if y_label else uv_absorbance_column_name)

            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.ylim(0, 2000)
            plt.xlim(6, 24)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # Adjust layout
        plt.tight_layout()

        # Save the plot for this group of files
        save_path = f"Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20250218_MSP1D1 sf GFP clean up{i // num_subplots + 1}.png"
        plt.savefig(save_path, dpi=1000)
        plt.show()
        plt.close()  # Close the figure to avoid displaying it

# Replace the list with the actual paths to your CSV files
file_paths = [
     "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20250206_MSP clean up\\20250213_OF-(+)MSP1D1 sf GFP Superdex200_1.csv"
   # "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\E3_mCherry Data\\20250214_OF_(+)MSP1E3D1 mCherry POPC 130_1 Rep2.csv",
  # "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\E3_mCherry Data\\20250214_OF_(+)MSP1E3D1 mCherry POPC 130_1 Rep3.csv"
    #"Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20241220\\20241220_OF_Alcohol dehydrogenase.csv",
   # "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20241220\\20241220_OF_Catalase.csv"
    #"Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20241213\\20241213_OF_MSP1D1_mCherry DMPC 90_1.csv",
   # "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20241213\\20241213_OF_MSP1D1_mCherry DMPC 100_1.csv"

]

# Provide custom legend names corresponding to each file path
legend_names = [
    '(+)MSP1D1 sf GFP clean up',
   # '(+)MSP1E3D1 POPC mCherry 2',
    #
   #'Alcohol Dehydrogenase',
    #'Catalase',
   # 'mCherry MSP1D1 DMPC nano 100:1'

]

# Specify custom axis labels if desired
custom_x_label = 'Retention Volume (ml)'
custom_y_label = 'A280 nm (mAU)'

plot_chromatograms(file_paths, legend_names, x_label=custom_x_label, y_label=custom_y_label)
