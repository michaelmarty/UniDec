import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file path for the Raw TIC
file_path = r'C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_carbonicanhydrase_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-21-09-53-49_unidecfiles\TIC_chrom.txt'

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

# Define color for Raw TIC
raw_color = '#001f4d'  # Dark blue for raw TIC

# --- Raw TIC only plot ---
plt.figure(figsize=(16, 9))

# Plot only the Raw TIC
data_raw = np.loadtxt(file_path)
time_raw = data_raw[:, 0]
intensity_raw = data_raw[:, 1]
plt.plot(time_raw, intensity_raw, '-', color=raw_color, lw=4, label='Raw TIC')

# Set labels and title
plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Intensity (# of ions)', labelpad=18)
plt.title('Raw Total Ion Chromatogram', pad=18)

# Format axes
ax_raw = plt.gca()
ax_raw.spines['top'].set_visible(False)
ax_raw.spines['right'].set_visible(False)
ax_raw.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_raw.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_raw.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])

# Set x-axis limits from 0 to 122 minutes
ax_raw.set_xlim(0, 122)

# Compute y-limits for Raw TIC
ymin_raw = np.min(intensity_raw)
ymax_raw = np.max(intensity_raw)

# Ensure ymin is non-negative
if ymin_raw < 0:
    ymin_raw = 0.0

# Set y-axis limits
ax_raw.set_ylim(ymin_raw, 60)

# Ensure 60 is included as a y-tick
yticks_raw = list(ax_raw.get_yticks())
if 60 not in yticks_raw:
    yticks_raw.append(60)
    yticks_raw = sorted(set(yticks_raw))
    ax_raw.set_yticks(yticks_raw)

# Save the Raw TIC plot as PNG and PDF
output_raw_png = os.path.join(os.path.dirname(file_path), 'raw_tic.png')
output_raw_pdf = os.path.join(os.path.dirname(file_path), 'raw_tic.pdf')
plt.tight_layout()
plt.savefig(output_raw_png, dpi=300, transparent=True)
plt.savefig(output_raw_pdf, dpi=300, transparent=True)
plt.close()

