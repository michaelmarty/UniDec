import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file paths for the three Decon TICs
file_paths = [
    r'C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_proteinmix_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-22-04-01-27_unidecfiles\TIC_decon.txt',
    r'C:\Users\bht442\Desktop\20251021\20251022_BHT_0.1mgmL_proteinmix_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-22-06-52-16_unidecfiles\TIC_decon.txt',
    r'C:\Users\bht442\Desktop\20251021\20251022_BHT_0.1mgmL_proteinmix_6inj_0.1gas_50inj_zoom_CDMS_C4_2_2025-10-22-08-56-16_unidecfiles\TIC_decon.txt'
]

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

# Define color scheme for the three datasets
colors = ['#001f4d', '#6A4C93', '#FF5733']  # Dark blue, muted purple, orange
labels = ['Decon TIC 1', 'Decon TIC 2', 'Decon TIC 3']

# Create a single figure for all datasets
plt.figure(figsize=(16, 9))

# Iterate over each file and plot the data
for i, file_path in enumerate(file_paths):
    # Load the data
    data = np.loadtxt(file_path)

    # Assuming the data has two columns: time and intensity
    time = data[:, 0]
    intensity = data[:, 1]

    # Plot the data
    plt.plot(time, intensity, '-', color=colors[i], lw=4, label=labels[i])

# Set labels and title
plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Intensity (# of ions)', labelpad=18)
plt.title('Total Ion Chromatogram Comparison', pad=18)

# Format axes
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
ax.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])

# Set x-axis to start at 0
ax.set_xlim(5, 15)

# Set y-axis limits to 0 to 140
ax.set_ylim(0, 140)

# Force the top of the y-axis to 140 and ensure it is a tick
top = 140.0
yticks = list(ax.get_yticks())
if top not in yticks:
    yticks.append(top)
    yticks = sorted(set(yticks))
    ax.set_yticks(yticks)

# Add legend without a surrounding box
ax.legend(frameon=False)

# Adjust layout to prevent cropping of labels
plt.tight_layout()

# Save the plot as a 300 DPI PNG and PDF with a transparent background in the same folder
output_path_png = os.path.join(os.path.dirname(file_paths[0]), 'decon_zoomed_tic_5_15.png')
output_path_pdf = os.path.join(os.path.dirname(file_paths[0]), 'decon_zoomed_tic_5_15.pdf')
plt.savefig(output_path_png, dpi=300, transparent=True)
plt.savefig(output_path_pdf, dpi=300, transparent=True)

plt.close()
