import os
import matplotlib.pyplot as plt
import numpy as np

# Define the file paths
file_paths = [
    r'Z:\Group Share\BHT\Agilent HPLC Data\20240730\20240730_BHT_GROEL_HADAMARD5BIT_5UL\20240730_BHT_GROEL_HADAMARD5BIT_5UL_VIAL100_TRACE.CSV',
    r'Z:\Group Share\BHT\Agilent HPLC Data\20240730\20240730_BHT_GROEL_HADAMARD5BIT_5UL\20240730_BHT_GROEL_HADAMARD5BIT_5UL_VIAL91_TRACE.CSV'
]

# Set font properties
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18

# Create a figure
plt.figure(figsize=(10, 6))

# Iterate over each file and plot the data
for file_path in file_paths:
    # Load the data with UTF-16 LE encoding and tab delimiter
    data = np.loadtxt(file_path, delimiter='\t', encoding='utf-16')

    # Assuming the data has two columns: time and absorbance
    time = data[:, 0]
    absorbance = data[:, 1]

    # Determine the label based on the file path
    label = "Vial 100" if "VIAL100" in file_path.upper() else "Vial 91"

    # Plot the data
    plt.plot(time, absorbance, label=label)

# Set labels
plt.xlabel('Time (minutes)')
plt.ylabel('Absorbance at 280 nm (mAU)')

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add a legend at the top left position with font size 11
plt.legend(loc='upper left', prop={'size': 11})

# Save the plot as a 300 DPI PNG with a transparent background in the same folder
output_path = os.path.join(os.path.dirname(file_paths[0]), 'chromatogram_comparison.png')
plt.savefig(output_path, dpi=300, transparent=True)

# Show the plot
plt.show()