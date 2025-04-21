import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.sans-serif'] = "Arial"




# Function to read CSV files and plot chromatograms with custom legend names and axis labels
def plot_chromatograms(file_paths, legend_names, x_label=None, y_label=None):

    plt.figure(figsize=(10, 6))

    for idx, file_path in enumerate(file_paths):
        # Read the CSV file, specifying the row containing column names (0-indexed)
        data = pd.read_csv(file_path, header=1)

        # Update these variable names based on the actual column names in your CSV
        volume_column_name = data.columns[0]  # The first column name
        uv_absorbance_column_name = data.columns[1]  # The second column name

        # Assuming the CSV file has the correct columns
        volume = data[volume_column_name]
        uv_absorbance = data[uv_absorbance_column_name]

        # Use the provided legend name for each plot
        label = legend_names[idx]

        # Plot chromatogram with a different colored line for each file
        # Bold the 6th file plot line
        if idx == 5:  # Remember that the index is 0-based, so 5 is the 6th file
            plt.plot(volume, uv_absorbance, label=label, linewidth=2)
        else:
            plt.plot(volume, uv_absorbance, label=label, linewidth=2)

    # Set the y-axis limit from 0 to 500
    plt.ylim(0, 500)

    plt.xlim(6, 24)

    # Add titles and labels
    plt.title('Fluorescent Nanodiscs')
    plt.xlabel(x_label if x_label else volume_column_name)
    plt.ylabel(y_label if y_label else uv_absorbance_column_name)
    # Add legend to the right outside of the plot area
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.legend(loc="upper right")

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)



    # Adjust layout
    plt.tight_layout()

    #plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20241023_Fluorescent naodiscs.png", dpi=1000)
    # Show the plot
    plt.show()


# Replace the list with the actual paths to your CSV files
file_paths = [
    "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\2024-10-03_OF_POPC control 65_1.csv",
    "Z:\\Group Share\\Seun\\LX-MS\\Biorad data\\20241022_OF_fnano 40_1.csv"

    # Add your file paths here
]

# Provide custom legend names corresponding to each file path
legend_names = [
    'Normal nanodisc control 65:1',
    'fnanodisc 40:1'

    # Add your custom names here
]

# Specify custom axis labels if desired
custom_x_label = 'Retention Volume (ml)'
custom_y_label = 'UV Absorbance (mAU)'

plot_chromatograms(file_paths, legend_names, x_label=custom_x_label, y_label=custom_y_label)

