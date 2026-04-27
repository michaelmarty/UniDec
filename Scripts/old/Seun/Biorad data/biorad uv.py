import pandas as pd
import matplotlib.pyplot as plt
import os


# Function to read CSV files and plot chromatograms
def plot_chromatograms(file_paths):
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
        # Label each plot with a number
        label = f'Chromatogram {idx + 1}'
        # Plot chromatogram with a different colored line for each file
        # Bold the 6th file plot line
        if idx == 5:  # Remember that the index is 0-based, so 5 is the 6th file
            plt.plot(volume, uv_absorbance, label=label, linewidth=2.5)
        else:
            plt.plot(volume, uv_absorbance, label=label, linewidth=1)
    # Set the y-axis limit from 0 to 2000
    plt.ylim(0, 600)
    # Add titles and labels
    plt.title('Chromatograms')
    plt.xlabel(volume_column_name)
    plt.ylabel(uv_absorbance_column_name)
    # Add legend to the right outside of the plot area
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # Show the plot
    plt.show()


# Replace the list with the actual paths to your 14 CSV files
file_paths = [
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\Biorad data\\20241022_OF_fnano 40_1.csv",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\Biorad data\\20241009_OF_fnano 50_1.csv",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\Biorad data\\2024-09-26_OF_fMSP1D1sfGFP_POPC.csv",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\Biorad data\\2024-10-03_OF_POPC fnano 80_1 whole mass.csv",

]

plot_chromatograms(file_paths)