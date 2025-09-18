import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Define the folder path
folder_path = r'Z:\Group Share\BHT\Agilent HPLC Data\20240801\20240801_BHT_GROEL_POSITIONTESTSEAT_7'

# Find all CSV files in the folder with "trace" in their filename
file_pattern = os.path.join(folder_path, '*trace*.CSV')
file_paths = glob.glob(file_pattern)

# Set font properties
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18

# Create a single figure
plt.figure(figsize=(10, 6))

# Initialize the absorbance offset
absorbance_offset = 0

# Iterate over each file and plot the data
for file_path in file_paths:
    # Load the data with UTF-16 LE encoding and tab delimiter
    data = np.loadtxt(file_path, delimiter='\t', encoding='utf-16')

    # Assuming the data has two columns: time and absorbance
    time = data[:, 0]
    absorbance = data[:, 1]

    # Extract the vial number from the filename
    vial_number = os.path.basename(file_path).split('_')[-2]

    # Add the absorbance offset to the absorbance values
    absorbance = absorbance + absorbance_offset

    # Plot the data
    plt.plot(time, absorbance, label=vial_number)

    # Increment the absorbance offset for the next dataset
    absorbance_offset += 10

# Set labels
plt.xlabel('Time (minutes)')
plt.ylabel('Absorbance at 280 nm (mAU)')

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limit to accommodate the staggered plots
ax.set_ylim(0, absorbance_offset + 10)
ax.set_xlim(0, 16)

# Add a legend
plt.legend(loc='upper right', prop={'size': 11})

# Save the plot as a 300 DPI PNG with a transparent background in the same folder
output_path = os.path.join(folder_path, 'chromatogram_comparison.png')
plt.savefig(output_path, dpi=300, transparent=True)

# Show the plot
plt.show()