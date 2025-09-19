import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file path
file_path = r'Z:\Group Share\BHT\Q Exactive HF Data\Walti HSP60\20250728\20250728_BHT_1mgmL_cHSP60_native_afterweekend_biospin_1_unidecfiles\20250728_BHT_1mgmL_cHSP60_native_afterweekend_biospin_1_rawdata.txt'
# Set font properties
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Create a new figure
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

# Make ticks bigger
ax.tick_params(axis='both', which='major', length=6, width=1.5)
ax.tick_params(axis='both', which='minor', length=6, width=1.5)

# Make the outside box lines (spines) thicker
for spine in ['left', 'bottom']:
    ax.spines[spine].set_linewidth(1.5)

# Set y-axis limits to 1.0
ax.set_ylim(0, 1.0)

# Set x-axis limits from 1000 to 16000
ax.set_xlim(1000, 16000)

# Set x-axis to display numbers without scientific notation
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Save the plot as a 300 DPI PNG with a transparent background in the same folder
output_path = os.path.join(os.path.dirname(file_path), 'mass_spectrum.png')
plt.savefig(output_path, dpi=300, transparent=True)

# Show the plot
plt.show()