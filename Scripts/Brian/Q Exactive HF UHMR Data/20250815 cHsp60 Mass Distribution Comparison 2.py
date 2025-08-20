import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file paths for all datasets to overlay
file_paths = [
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250814\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_1_2025-08-14-02-39-41_unidecfiles\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_1_2025-08-14-02-39-41_mass.txt',
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250814\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_2_2025-08-14-03-38-00_unidecfiles\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_2_2025-08-14-03-38-00_mass.txt',
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250814\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_3_2025-08-14-04-21-16_unidecfiles\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_3_2025-08-14-04-21-16_mass.txt',
    r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250814\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_4_2025-08-14-04-55-51_unidecfiles\20250814_BHT_cHsp60_0.05mgmL_directinfusion_CMDS_4_2025-08-14-04-55-51_mass.txt'
]

# Labels and colors for each dataset
labels = ['Replicate 1', 'Replicate 2', 'Replicate 3', 'Replicate 4']
colors = ['black', 'red', 'blue', 'green']

# Set font and figure size for the plot
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28
plt.figure(figsize=(16, 9))

# Loop through each file and plot its data
for i, file_path in enumerate(file_paths):
    # Load data from file
    data = np.loadtxt(file_path)
    mass = data[:, 0] / 1000  # Convert mass to kDa
    intensity = data[:, 1]
    intensity = intensity / np.max(intensity)  # Normalize intensity

    # Plot the data with label and color
    plt.plot(mass, intensity, label=labels[i], color=colors[i])

# Set axis labels
plt.xlabel('Mass (kDa)')
plt.ylabel('Normalized Intensity')

# Get current axis for further formatting
ax = plt.gca()

# Hide top and right spines (box lines)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis and x-axis limits
ax.set_ylim(0, 1.0)
ax.set_xlim(0, 1200)

# Format x-axis to show plain numbers (not scientific notation)
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Show legend for the datasets
plt.legend()

# Save the figure to the same directory as the first file
output_path = os.path.join(os.path.dirname(file_paths[0]), 'mass_distribution_comparison.png')
plt.savefig(output_path, dpi=300, transparent=True)

# Display the plot
plt.show()