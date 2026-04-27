import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import numpy as np

rcParams['font.sans-serif'] = "Arial"
# === FILES TO PLOT ===
files = [
    r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251029\20251029_EC_D1sfCherry_short_1to80_DMPC_BCA_CDMS_2025-10-29-10-36-35_unidecfiles\20251029_EC_D1sfCherry_short_1to80_DMPC_BCA_CDMS_2025-10-29-10-36-35_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251029\20251029_EC_D1sfGFP_1to80_DMPC_BCA_CDMS_2025-10-29-09-46-19_unidecfiles\20251029_EC_D1sfGFP_1to80_DMPC_BCA_CDMS_2025-10-29-09-46-19_mass.txt"
    # Add your file paths here
]

# === LEGEND LABELS ===
labels = [
"MSP1D1 1:81 DMPC",
    "MSP1D1-sfCherry 1:79  DMPC",
"MSP1D1-sfGFP 1:88 DMPC",
    # Add your legend labels here
]

# === FIGURE SETUP ===
fig, axes = plt.subplots(len(files), 1, sharex=True, figsize=(6, 4), dpi=300)
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

# === LOOP THROUGH FILES ===
for i, (file, label) in enumerate(zip(files, labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize to 100%
    df["Mass_kDa"] = df["Mass_Da"] / 1000
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    ax = axes[i]
    offset = i * y_offset_step
    color = colors[i % len(colors)]

    # Plot data
    ax.plot(df["Mass_kDa"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)

    # Styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(100, 300)
    ax.set_ylim(offset, offset + 100)
    ax.set_yticks([offset, offset + 50, offset + 100])
    ax.set_yticklabels(y_labels)

    # Add legend to the left **inside the axes**, but offset from y-axis
    ax.legend(loc="upper left", bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=9)

    # Only show x-axis labels on the bottom plot
    if i < len(files) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Mass (kDa)")

# Common y-axis label (moved slightly left)
fig.text(0.03, 0.5, "Relative Intensity (%)", va="center", rotation="vertical")

# # Transparent background
# fig.patch.set_alpha(0.0)
# for ax in axes:
#     ax.set_facecolor("none")

# === FIX CLIPPING & SPACING ===
# plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.4)
path=r"D:\figs"
plt.savefig(os.path.join(path, "fig1try100.pdf"))
plt.show()
