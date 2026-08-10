import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file paths
file_paths = [
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250327\Control CDMS\20250327_BHT_HSP60_native_1_unidecfiles\20250327_BHT_HSP60_native_1_mass.txt'
]

# Set font properties
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Create a single figure
plt.figure(figsize=(16, 9))

# Iterate over each file and plot the data
for i, file_path in enumerate(file_paths):
    # Load the data
    data = np.loadtxt(file_path)

    # Assuming the data has two columns: mass and intensity
    mass = data[:, 0] / 1000  # Convert mass to kDa
    intensity = data[:, 1]

    # Normalize the intensity
    intensity = intensity / np.max(intensity)

    # Set the label based on the file index
    label = 'Nanoflow' if i == 0 else 'Static' if i == 1 else 'New Data'

    # Plot the data with black lines
    plt.plot(mass, intensity, label=label, color='black')

# Set labels
plt.xlabel('Mass (kDa)')
plt.ylabel('Normalized Intensity')

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limits to end at 1.0
ax.set_ylim(0, 1.0)

# Set x-axis limits to end at 1,000 kDa
ax.set_xlim(0, 1000)

# Set x-axis to display numbers without scientific notation
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Save the plot as a 300 DPI PNG with a transparent background in the same folder
output_path = os.path.join(os.path.dirname(file_paths[0]), 'mass_distribution_comparison.png')
plt.savefig(output_path, dpi=300, transparent=True)

# Show the plot
plt.show()