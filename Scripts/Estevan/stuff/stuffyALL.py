import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
x_range = (100, 300)  # x-axis limits

# === FILES TO PLOT ===
files = [
r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251016\20251016_EC_(-)D1T0-sfcherryLong_DMPC150x_CDMS_2025-10-16-11-46-29_unidecfiles\20251016_EC_(-)D1T0-sfcherryLong_DMPC150x_CDMS_2025-10-16-11-46-29_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251203\20251203EC_D1sfGFP_1to80_DMPCcdms_2025-12-03-11-33-15_unidecfiles\20251203EC_D1sfGFP_1to80_DMPCcdms_2025-12-03-11-33-15_mass.txt"
    # Add your file paths here
]

# === MAIN LEGEND LABELS (one per file) ===
labels = [
"MSP1D1",
    "MSP1D1-sfCherryS",
"MSP1D1-sfGFP",
    # Add your legend labels here
]

# === SECOND LEGEND LABELS (one per file) ===
second_labels = [
"162 ± 4 DMPC per Nanodisc",
    "166 ± 10 DMPC per Nanodisc",
    "162 ± 7 DMPC per Nanodisc",
    # Add your second legend labels here, e.g. "FWHM", "Peak"
]

# === FIGURE SETUP ===
fig, axes = plt.subplots(len(files), 1, sharex=True, figsize=(9, 5), dpi=300)
if len(files) == 1:
    axes = [axes]

# === CUSTOM COLORS (dark viridis subset) ===
colors = [
    "#0F0C08",  # deep purple
    "#CC79A7",  # dark blue
    "#009E73"   # blue-green
]

# === Y-OFFSET SETTINGS ===
y_offset_step = 120
y_labels = ["0", "%", "100"]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# --- LOOP THROUGH FILES ---
for i, (file, label, second_label) in enumerate(zip(files, labels, second_labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    ax = axes[i]
    offset = i * y_offset_step
    color = colors[i % len(colors)]

    # Plot the main curve
    line, = ax.plot(df["Mass_kDa"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)

    # --- Find peak for second legend ---
    peak_idx = df["Percent"].idxmax()
    peak_mass = df["Mass_kDa"].iloc[peak_idx]
    peak_value = df["Percent"].iloc[peak_idx] + offset

    # Styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(x_range)  # configurable x-axis range
    ax.set_ylim(offset, offset + 100)
    ax.set_yticks([offset, offset + 50, offset + 100])
    ax.set_yticklabels(y_labels)

    # --- First legend (curve name, upper left) ---
    curve_legend = ax.legend([line], [label], loc="upper left",
                             frameon=False, fontsize=font_size)
    ax.add_artist(curve_legend)  # explicitly add first legend

    # --- Second legend (near peak) ---
    # Convert peak position to axes fraction for proper placement
    x_frac = (peak_mass - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])
    y_frac = (peak_value - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    second_handle = Line2D([0], [0], color=color, lw=1.5)
    second_legend = ax.legend([second_handle], [second_label],
                              loc='center', bbox_to_anchor=(x_frac + 0.05, y_frac),
                              frameon=False, fontsize=font_size)
    ax.add_artist(second_legend)

    # Only show x-axis labels on the bottom plot
    if i < len(files) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Mass (kDa)", fontname=font_name, fontsize=font_size)

# Common y-axis label
fig.text(0.05, 0.5, "Relative Intensity (%)", va="center", rotation="vertical",
         fontname=font_name, fontsize=font_size)

# Adjust figure spacing
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.6)

plt.show()
