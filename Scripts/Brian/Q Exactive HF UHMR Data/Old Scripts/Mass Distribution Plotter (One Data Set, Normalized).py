import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Prompt the user for the file path and strip quotes
file_path = input("Enter the file path (with quotes): ").strip('"').strip("'")

# Set font properties and styling
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28
plt.rcParams['xtick.major.size'] = 12
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 12
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
plt.rcParams['axes.labelpad'] = 18
plt.rcParams['axes.titlepad'] = 18
plt.rcParams['axes.linewidth'] = plt.rcParams['xtick.major.width']
plt.rcParams['pdf.fonttype'] = 42 # Ensure fonts are saved as editable text in PDFs

# Create a new figure
plt.figure(figsize=(16, 9))

# Load the data
data = np.loadtxt(file_path)

# Assuming the data has two columns: mass and intensity
mass = data[:, 0]
intensity = data[:, 1]

# Convert mass from Da to kDa
mass = mass / 1000

# Normalize the intensity
intensity = intensity / np.max(intensity)

# Plot the data with black lines
plt.plot(mass, intensity, color='black')

# Set labels
plt.xlabel(r'Mass (kDa)')
plt.ylabel('Normalized Intensity')

# Remove the top and right lines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limits to 1.0
ax.set_ylim(0, 1.0)

# Set x-axis limits based on the dataset
ax.set_xlim(np.min(mass), np.max(mass))

# Set x-axis to display numbers without scientific notation
formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# Automatically adjust the layout to prevent label cropping
plt.tight_layout()

# Save the plot as a 300 DPI PNG with a transparent background in the same folder
output_path_png = os.path.join(os.path.dirname(file_path), 'mass_distribution.png')
plt.savefig(output_path_png, dpi=300, transparent=True)

# Save the plot as a PDF in the same folder
output_path_pdf = os.path.join(os.path.dirname(file_path), 'mass_distribution.pdf')
plt.savefig(output_path_pdf, dpi=300)

# Print the output paths
print(f"Plot saved as PNG at: {output_path_png}")
print(f"Plot saved as PDF at: {output_path_pdf}")

# Show the plot
plt.show()