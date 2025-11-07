import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file path for the Raw TIC
file_path = r'C:\Users\bht442\Desktop\20251030\20251030_BHT_POSB_10XDiluted_HEKcelllysate_1INJ_C4_1_unidecfiles\Raw TIC.txt'

# Define file paths for the additional replicates
file_path_2 = r'C:\Users\bht442\Desktop\20251030\20251030_BHT_POSB_10XDiluted_HEKcelllysate_1INJ_C4_2_unidecfiles\Raw TIC.txt'
file_path_3 = r'C:\Users\bht442\Desktop\20251030\20251030_BHT_POSB_10XDiluted_HEKcelllysate_1INJ_C4_3_unidecfiles\Raw TIC.txt'

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

# Define color for Raw TIC
raw_color = '#001f4d'  # Dark blue for raw TIC
# Define colors for the additional replicates
replicate_2_color = '#004d00'  # Dark green for replicate 2
replicate_3_color = '#660000'  # Dark red for replicate 3

# --- Raw TIC only plot ---
plt.figure(figsize=(16, 9))

# Plot only the Raw TIC
data_raw = np.loadtxt(file_path)
time_raw = data_raw[:, 0]
intensity_raw = data_raw[:, 1]

# Normalize the intensity values
intensity_raw = (intensity_raw - np.min(intensity_raw)) / (np.max(intensity_raw) - np.min(intensity_raw))

plt.plot(time_raw, intensity_raw, '-', color=raw_color, lw=4, label='Raw TIC')

# --- Overlay plot for Raw TIC and replicates ---
plt.figure(figsize=(16, 9))

# Plot the first Raw TIC
data_raw = np.loadtxt(file_path)
time_raw = data_raw[:, 0]
intensity_raw = data_raw[:, 1]
intensity_raw = (intensity_raw - np.min(intensity_raw)) / (np.max(intensity_raw) - np.min(intensity_raw))
plt.plot(time_raw, intensity_raw, '-', color=raw_color, lw=4, label='Raw TIC (Replicate 1)')

# Plot the second replicate
data_replicate_2 = np.loadtxt(file_path_2)
time_replicate_2 = data_replicate_2[:, 0]
intensity_replicate_2 = data_replicate_2[:, 1]
intensity_replicate_2 = (intensity_replicate_2 - np.min(intensity_replicate_2)) / (np.max(intensity_replicate_2) - np.min(intensity_replicate_2))
plt.plot(time_replicate_2, intensity_replicate_2, '-', color=replicate_2_color, lw=4, label='Raw TIC (Replicate 2)')

# Plot the third replicate
data_replicate_3 = np.loadtxt(file_path_3)
time_replicate_3 = data_replicate_3[:, 0]
intensity_replicate_3 = data_replicate_3[:, 1]
intensity_replicate_3 = (intensity_replicate_3 - np.min(intensity_replicate_3)) / (np.max(intensity_replicate_3) - np.min(intensity_replicate_3))
plt.plot(time_replicate_3, intensity_replicate_3, '-', color=replicate_3_color, lw=4, label='Raw TIC (Replicate 3)')

# Set labels and title
plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Normalized Intensity', labelpad=18)
plt.title('Overlay of Raw Total Ion Chromatograms', pad=18)

# Add legend
plt.legend(loc='upper left', fontsize=18)  # Changed legend location to 'upper left'

# Format axes
ax_overlay = plt.gca()
ax_overlay.spines['top'].set_visible(False)
ax_overlay.spines['right'].set_visible(False)
ax_overlay.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_overlay.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_overlay.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])
plt.xlim(left=0, right=20)  # Ensure x-axis starts at 0 and ends at 20
plt.ylim(bottom=0, top=1)  # Ensure y-axis starts at 0 and ends at 1

# Save the overlay plot as PNG and PDF
output_overlay_png = os.path.join(os.path.dirname(file_path), 'overlay_tic.png')
output_overlay_pdf = os.path.join(os.path.dirname(file_path), 'overlay_tic.pdf')
plt.tight_layout()
plt.savefig(output_overlay_png, dpi=300, transparent=True)
plt.savefig(output_overlay_pdf, dpi=300, transparent=True)
plt.close()
