import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import os

# Define the file path for the Raw TIC
file_path = r'C:\Users\bht442\Desktop\20251030\20251030_BHT_POSB_10XDiluted_HEKcelllysate_1INJ_C4_1_unidecfiles\Raw TIC.txt'

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

# --- Raw TIC only plot ---
plt.figure(figsize=(16, 9))

# Plot only the Raw TIC
data_raw = np.loadtxt(file_path)
time_raw = data_raw[:, 0]
intensity_raw = data_raw[:, 1]

# Normalize the intensity values
intensity_raw = (intensity_raw - np.min(intensity_raw)) / (np.max(intensity_raw) - np.min(intensity_raw))

plt.plot(time_raw, intensity_raw, '-', color=raw_color, lw=4, label='Raw TIC')

# Set labels and title
plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Normalized Intensity', labelpad=18)
plt.title('Raw Total Ion Chromatogram', pad=18)

# Format axes
ax_raw = plt.gca()
ax_raw.spines['top'].set_visible(False)
ax_raw.spines['right'].set_visible(False)
ax_raw.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_raw.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
ax_raw.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])
plt.xlim(left=0, right=20)  # Ensure x-axis starts at 0 and ends at 20
plt.ylim(bottom=0, top=1)  # Ensure y-axis starts at 0 and ends at 1

# Save the Raw TIC plot as PNG and PDF
output_raw_png = os.path.join(os.path.dirname(file_path), 'raw_tic.png')
output_raw_pdf = os.path.join(os.path.dirname(file_path), 'raw_tic.pdf')
plt.tight_layout()
plt.savefig(output_raw_png, dpi=300, transparent=True)
plt.savefig(output_raw_pdf, dpi=300, transparent=True)
plt.close()

