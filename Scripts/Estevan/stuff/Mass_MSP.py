import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
x_range = (0, 500)  # x-axis limits

# === FILES TO PLOT ===
files = [
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfCherry - done\Long - done\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_unidecfiles\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mCherry - done\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_unidecfiles\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mBlueberry - done\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_unidecfiles\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\Halo-Tag - done\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_unidecfiles\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfGFP - done\E3T0\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_unidecfiles\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_mass.txt"
]

# === MAIN LEGEND LABELS (one per file) ===
labels = [
    "MSP1D1-sfCherryL",
    "MSP1D1-mCherry",
    "MSP1D1-mBlueberry",
    "MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP"
]

# === SECOND LEGEND LABELS (one per file) ===
second_labels = [
    "185 ± 16 DMPC per Nanodisc",
    "154 ± 9 DMPC per Nanodisc",
    "156 ± 7 DMPC per Nanodisc",
    "146 ± 14 DMPC per Nanodisc",
    "360 ± 20 DMPC per Nanodisc"
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(9, 5), dpi=300)

# === CUSTOM COLORS (dark viridis subset) ===
colors = [
    "#CC79A7",  # deep purple
    "#CC79A7",  # dark blue
    "#2C7FB8",
    "#7A7A7A",  # blue-green
    "#00A9A5"
]

# === Y-OFFSET SETTINGS ===
y_offset_step = 120
y_labels = ["0", "%", "100"]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# For legend handles
handles = []

# Track global min/max for nicer axis limits if desired
min_y_display = np.inf
max_y_display = -np.inf

# --- LOOP THROUGH FILES AND PLOT ON THE SAME AXIS ---
for i, (file, label, second_label) in enumerate(zip(files, labels, second_labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize to percent
    df["Mass_kDa"] = df["Mass_Da"] / 1000
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    offset = i * y_offset_step
    color = colors[i % len(colors)]

    # Plot the main curve on the single axis
    line, = ax.plot(df["Mass_kDa"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)
    handles.append(line)

    # Find peak for placing the second label as an annotation
    peak_idx = df["Percent"].idxmax()
    peak_mass = df["Mass_kDa"].iloc[peak_idx]
    peak_value = df["Percent"].iloc[peak_idx] + offset

    # Annotate the second label near the peak
    # Slight horizontal offset to avoid overlapping the line
    xlim = ax.get_xlim() if ax.has_data() else x_range
    # compute fractional offset in axis coordinates if needed; simpler: small shift in kDa
    x_shift = (xlim[1] - xlim[0]) * 0.02
    ax.text(peak_mass + x_shift, peak_value, second_label,
            fontsize=font_size * 0.9, color=color, va='center', ha='left', alpha=0.9)

    # Update min/max display extents for y axis
    min_y_display = min(min_y_display, offset)
    max_y_display = max(max_y_display, offset + 101.5)  # each trace displays 0..100 (+ small margin)

# Single legend for all traces
# To avoid a giant legend when many traces exist, limit entries if needed (your current list is small)
ax.legend(handles, labels, loc='upper right', frameon=False, fontsize=font_size)

# Axis styling and limits
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)  # configurable x-axis range

# Set y-limits to encompass all stacked traces
if min_y_display == np.inf:
    min_y_display = 0
if max_y_display == -np.inf:
    max_y_display = y_offset_step * (len(files) - 1) + 101.5
ax.set_ylim(min_y_display - 5, max_y_display + 5)

# Create y-ticks showing 0, 50, 100 for each offset level
yticks = []
yticklabels = []
for i in range(len(files)):
    offset = i * y_offset_step
    yticks.extend([offset, offset + 50, offset + 100])
    # For clarity we show only the repeating labels "0", "%", "100" for each trace row
    yticklabels.extend(y_labels)

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

# Labels
ax.set_xlabel("Mass (kDa)", fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)

plt.show()