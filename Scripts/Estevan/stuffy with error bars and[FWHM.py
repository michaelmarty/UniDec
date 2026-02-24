import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === FILES TO PLOT ===
files = [
r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251106\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_CDMS_2025-11-06-11-40-25_unidecfiles\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_CDMS_2025-11-06-11-40-25_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251105\20251105_EC_D1sfGFP1to80_DMPC_CDMSlowerspray_2025-11-05-11-20-47_unidecfiles\20251105_EC_D1sfGFP1to80_DMPC_CDMSlowerspray_2025-11-05-11-20-47_mass.txt",
r"Z:\Group Share\Estevan\Q Exactive Data\20251111\20251111_EC_D1mblueberry_1to80_DMPC_CDMS_2025-11-11-11-03-38_unidecfiles\20251111_EC_D1mblueberry_1to80_DMPC_CDMS_2025-11-11-11-03-38_mass.txt",
   # r"Z:\Group Share\Estevan\Q Exactive Data\20251112\20251212_EC_E3sfGFP1to_175_DMPC_cdms_2025-11-12-02-43-19_unidecfiles\20251212_EC_E3sfGFP1to_175_DMPC_cdms_2025-11-12-02-43-19_mass.txt",
    #r"Z:\Group Share\Estevan\Q Exactive Data\20251022\20251022_D1sfGFPminus_1to150_DMPC_CDMS_2025-10-22-02-34-17_unidecfiles\20251022_D1sfGFPminus_1to150_DMPC_CDMS_2025-10-22-02-34-17_mass.txt"
    # "your_file_1.txt",
    # "your_file_2.txt",
    # "your_file_3.txt"
]

# === LEGEND LABELS ===
labels = [
"MSP1D1",
    "MSP1D1-sfCherry",
"MSP1D1-sfGFP",
     "MSP1D1-mBlueberry2",
     #"MSP1E3D1-sfGFP",
    # "Sample C"
]

# === SECOND LEGEND LABELS (with ± sign) ===
second_labels = [
    "162 ± 4 DMPC per ND",
    "162 ± 3 DMPC per ND",
    "164 ± 6 DMPC per ND",
    "144 ± 6 DMPC per ND",
    #"378 ± 30 DMPC per ND",
    # or "FWHM \u00B1 5.2 kDa" if your editor can't handle ±
]

# === FIGURE SETUP ===
fig, axes = plt.subplots(len(files), 1, sharex=True, figsize=(6, 5), dpi=300)
if len(files) == 1:
    axes = [axes]

# === FONT CONFIG ===
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10
})

# === CUSTOM COLORS (dark viridis subset) ===
colors = [
    "#0F0C08",  # deep purple
    "#CC79A7",  # dark blue
    "#009E73",  # blue-green
    "#56B4E9",
    "#CC79A7",
]

# === Y-OFFSET SETTINGS ===
y_offset_step = 120
y_labels = ["0", "%", "100"]

# === LOOP THROUGH FILES ===
for i, (file, label, second_label) in enumerate(zip(files, labels, second_labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize to 100%
    df["Mass_kDa"] = df["Mass_Da"] / 1000
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    ax = axes[i]
    offset = i * y_offset_step
    color = colors[i % len(colors)]

    # Plot data
    line, = ax.plot(df["Mass_kDa"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)

    # Styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(100, 400)
    ax.set_ylim(offset, offset + 100)
    ax.set_yticks([offset, offset + 50, offset + 100])
    ax.set_yticklabels(y_labels)

    # === MAIN LEGEND (for curve name) ===
    ax.legend(loc="center left", bbox_to_anchor=(1.15, 0.5), frameon=False)

    # === SECOND "LEGEND" (text label for FWHM) ===
    # Find approximate x-position near the right side of the peak
    peak_x = df.loc[df["Percent"].idxmax(), "Mass_kDa"]
    right_x = peak_x + 70  # shift right
    ax.text(
        right_x, offset + 80,  # position
        second_label,
        color=color,
        fontsize=10,
        va="center",
        ha="left"
    )

    # Only show x-axis labels on the bottom plot
    if i < len(files) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Mass (kDa)")

# Common y-axis label
fig.text(0.03, 0.5, "Relative Intensity (%)", va="center", rotation="vertical")

# Transparent background
fig.patch.set_alpha(0.0)
for ax in axes:
    ax.set_facecolor("none")

# === FIX CLIPPING & SPACING ===
plt.subplots_adjust(left=0.2, right=0.88, bottom=0.15, top=0.95, hspace=0.45)
plt.tight_layout(rect=[0.05, 0, 0.9, 1])

plt.show()
