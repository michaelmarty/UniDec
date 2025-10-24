import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file paths
file_paths = [
    r'C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_carbonicanhydrase_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-21-09-53-49_unidecfiles\TIC_chrom.txt',
    r'C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_carbonicanhydrase_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-21-09-53-49_unidecfiles\TIC_decon.txt'
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

# Define color scheme
raw_color = '#001f4d'      # Dark blue for raw TIC
demux_color = '#6A4C93'    # Muted purple for demultiplexed TIC

# --- Zoomed-in plot (10 to 14 minutes) ---
plt.figure(figsize=(16, 9))

# Plot both Raw and Demultiplexed TIC within the zoom range
zoom_min, zoom_max = 10, 14
for i, file_path in enumerate(file_paths):
    data = np.loadtxt(file_path)
    time = data[:, 0]
    intensity = data[:, 1]

    # Apply zoom range
    mask = (time >= zoom_min) & (time <= zoom_max)
    time_zoom = time[mask]
    intensity_zoom = intensity[mask]

    # Choose color and label
    color = raw_color if i == 0 else demux_color
    label = 'Raw TIC' if i == 0 else 'Demultiplexed TIC'

    # Plot the zoomed data
    plt.plot(time_zoom, intensity_zoom, '-', color=color, lw=4, label=label)

# Set labels and title
plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Intensity (# of ions)', labelpad=18)
plt.title('Total Ion Chromatogram (Zoomed: 10-14 min)', pad=18)

# Format axes
ax_zoom = plt.gca()
ax_zoom.spines['top'].set_visible(False)
ax_zoom.spines['right'].set_visible(False)
ax_zoom.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_zoom.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_zoom.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])

# Set x-axis limits to the zoom range
ax_zoom.set_xlim(zoom_min, zoom_max)

# Compute y-limits for zoomed-in plot
ymin_zoom = min(np.min(data[:, 1][(data[:, 0] >= zoom_min) & (data[:, 0] <= zoom_max)]) for data in [np.loadtxt(file) for file in file_paths])
ymax_zoom = max(np.max(data[:, 1][(data[:, 0] >= zoom_min) & (data[:, 0] <= zoom_max)]) for data in [np.loadtxt(file) for file in file_paths])

# Set y-axis limits
ax_zoom.set_ylim(ymin_zoom, max(60, ymax_zoom))

# Ensure 60 is included as a y-tick
yticks = list(ax_zoom.get_yticks())
if 60 not in yticks:
    yticks.append(60)
    yticks = sorted(set(yticks))
    ax_zoom.set_yticks(yticks)

# Add legend
ax_zoom.legend(frameon=False)

# Save the zoomed-in plot as PNG and PDF
output_zoom_png = os.path.join(os.path.dirname(file_paths[0]), 'zoomed_tic_10_14.png')
output_zoom_pdf = os.path.join(os.path.dirname(file_paths[0]), 'zoomed_tic_10_14.pdf')
plt.tight_layout()
plt.savefig(output_zoom_png, dpi=300, transparent=True)
plt.savefig(output_zoom_pdf, dpi=300, transparent=True)
plt.close()

