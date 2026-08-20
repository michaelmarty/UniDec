import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file paths
file_paths = [
    r'C:\Users\Brian Tran\University of Arizona\martylab - Documents\Presentations\Brian\Candidacy Exam\Figures\20250327_BHT_HSP60_CDMS_long_2_scan1812_unidecfiles\20250327_BHT_HSP60_CDMS_long_2_scan1812_rawdata.txt',
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250327\Control CDMS\20250327_BHT_HSP60_native_1_unidecfiles\20250327_BHT_HSP60_native_1_rawdata.txt'
]

# Set font properties
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Iterate over each file and plot the data
for i, file_path in enumerate(file_paths):
    # Create a new figure for each file
    plt.figure(figsize=(16, 9))

    # Load the data
    data = np.loadtxt(file_path)

    # Assuming the data has two columns: m/z and intensity
    mz = data[:, 0]
    intensity = data[:, 1]

    # Normalize the intensity
    intensity = intensity / np.max(intensity)

    # Plot the data with black lines
    plt.plot(mz, intensity, color='black')

    # Set labels
    plt.xlabel(r'$m/z$')
    plt.ylabel('Normalized Intensity')

    # Remove the top and right lines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set y-axis limits to 1.0
    ax.set_ylim(0, 1.0)

    # Set x-axis limits from 1000 to 20000
    ax.set_xlim(1000, 20000)

    # Set x-axis to display numbers without scientific notation
    formatter = ScalarFormatter(useOffset=False, useMathText=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)

    # Save the plot as a 300 DPI PNG with a transparent background in the same folder
    output_path = os.path.join(os.path.dirname(file_paths[0]), f'mass_spectrum_{i}.png')
    plt.savefig(output_path, dpi=300, transparent=True)

# Show all plots simultaneously
plt.show()