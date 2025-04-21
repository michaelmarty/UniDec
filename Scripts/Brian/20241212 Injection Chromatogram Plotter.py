import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Define the folder path
folder_path = r'Z:\Group Share\BHT\Agilent HPLC Data\20240806'

# Find all CSV files in the folder with "trace" in their filename
file_pattern = os.path.join(folder_path, '*trace*.CSV')
file_paths = glob.glob(file_pattern)

# Iterate over each file and plot the data
for file_path in file_paths:
    # Load the data with UTF-16 LE encoding and tab delimiter
    data = np.loadtxt(file_path, delimiter='\t', encoding='utf-16')

    # Assuming the data has two columns: time and absorbance
    time = data[:, 0]
    absorbance = data[:, 1]

    # Set font properties
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 18

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(time, absorbance)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Absorbance at 280 nm (mAU)')

    # Remove the top and right lines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set y-axis limit to 100 mAU and x-axis limit to 16 minutes
    ax.set_ylim(0, 30)
    ax.set_xlim(0, 30)

    # Save the plot as a 300 DPI PNG with a transparent background in the same folder
    output_path = os.path.join(os.path.dirname(file_path), os.path.splitext(os.path.basename(file_path))[0] + '.png')
    plt.savefig(output_path, dpi=300, transparent=True)

# Show all plots at the same time
plt.show()