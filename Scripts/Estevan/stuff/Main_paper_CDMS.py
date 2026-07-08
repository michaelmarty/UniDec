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
r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Long - done\20251208_EC_D1cherryLong_1to90_DMPC_rep2_cdms_2025-12-08-12-09-32_unidecfiles\20251208_EC_D1cherryLong_1to90_DMPC_rep2_cdms_2025-12-08-12-09-32_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\mCherry - done\20251202_EC_MSP1D1mCherry_CDMS_rep1_unidecfiles\20251202_EC_MSP1D1mCherry_CDMS_rep1_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\mBlueberry - done\20251111_EC_MSP1D1mBlueberry__CMS_rep2_unidecfiles\20251111_EC_MSP1D1mBlueberry__CMS_rep2_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\Halo-Tag - done\20251204_EC_D1Halo_1to80_DMPC_rep2cdms_2025-12-04-09-52-12_unidecfiles\20251204_EC_D1Halo_1to80_DMPC_rep2cdms_2025-12-04-09-52-12_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfGFP - done\E3T0\20251205_EC_MSP1E3D1sfGFP_CDMS_rep1_unidecfiles\20251205_EC_MSP1E3D1sfGFP_CDMS_rep1_mass.txt",
r"Z:\Group Share\Estevan\Q Exactive Data\Paper Revision E3T3\20251205_EC_E3t3_1to150DMPC_cdms_2025-12-05-11-07-52_unidecfiles\20251205_EC_E3t3_1to150DMPC_cdms_2025-12-05-11-07-52_mass.txt"
    # Add your file paths here
]

# === MAIN LEGEND LABELS (one per file) ===
labels = [
"MSP1D1-sfCherryL",
    "MSP1D1-mCherry",
"MSP1D1-mBlueberry2",
"MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP",
    "MSP1E3D1"
    # Add your legend labels here
]

# === SECOND LEGEND LABELS (one per file) ===
second_labels = [
"166 ± 23 DMPC per Nanodisc",
    "159 ± 10 DMPC per Nanodisc",
    "180 ± 16 DMPC per Nanodisc",
    "146 ± 14 DMPC per Nanodisc",
    "360 ± 20 DMPC per Nanodisc",
"364 ± 47 DMPC per Nanodisc",
    # Add your second legend labels here, e.g. "FWHM", "Peak"
]

# === FIGURE SETUP ===
fig, axes = plt.subplots(len(files), 1, sharex=True, figsize=(6, 5), dpi=300)
if len(files) == 1:
    axes = [axes]

# === CUSTOM COLORS (dark viridis subset) ===
colors = [
    "#CC79A7",  # deep purple
    "#CC79A7",  # dark blue
    "#2C7FBB",
    "#7A7A7A",# blue-green
   "#00A9A5",
"#3A3A3A"
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
    ax.set_ylim(offset, offset + 101.5)
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
    second_handle = Line2D([0], [0], color=color, lw=0.0)
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
