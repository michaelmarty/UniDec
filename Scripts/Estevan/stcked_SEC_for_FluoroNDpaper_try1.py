import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

# -----------------------
# Settings
# -----------------------
font_size = 10
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = font_size

# === FILES TO PLOT ===
files = [
    r"D:\FLND AKTA files\XLSX chromatograms\D1_1to80_DMPC.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\D1sfGFP_1to80_DMPC.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\D1sfCherrySHORT_1to80_DMPC.xlsx",
    r"D:\FLND AKTA files\New folder\MSP1d1_CHERRYLONG.xlsx",
    r"D:\FLND AKTA files\New folder\mcherry.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\d1bLUEBERRY_1TO80_dmpc.xlsx",
    r"D:\FLND AKTA files\XLSX chromatograms\HaloTag.xlsx",
    r"D:\FLND AKTA files\New folder\egfp.xlsx",
]

# === LEGEND LABELS ===
sample_labels = [
    "MSP1D1",
    "MSP1D1-sfGFP",
    "MSP1D1-sfCherryshort",
    "MSP1D1-sfCherryLong",
    "MSP1D1-mCherry",
    "MSP1D1-mBlueberry2",
    "MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP"
]

# === CUSTOM HEX COLORS ===
colors = [
    "#333333",  # charcoal
    "#00A9A5",  # green
    "#C76BA0",  # pink
    "#C76BA0",  # blue
    "#C76BA0",  # grey
    "#2C7FB8",
    "#7A7A7A",
    "#00A9A5",
]

# -----------------------
# Read Excel starting at A4/B4
# -----------------------
def read_sec_file(file_path):
    df = pd.read_excel(file_path, skiprows=3, header=None, usecols="A:B")
    df.columns = ["ml", "mAU"]
    df = df.dropna()
    return df

# -----------------------
# Plot setup
# -----------------------
fig, axes = plt.subplots(len(files), 1, figsize=(8, 10), sharex=True, dpi=300)
if len(files) == 1:
    axes = [axes]

for i, (file, label) in enumerate(zip(files, sample_labels)):
    df = read_sec_file(file)
    df["mAU"] = df["mAU"] / df["mAU"].max()  # normalize 0–1

    axes[i].plot(df["ml"], df["mAU"], color=colors[i % len(colors)], linewidth=1.8, label=label)

    # Axis formatting
    axes[i].set_ylim(0, 1.05)
    axes[i].set_xlim(10, 20)  # cutoff at 25 mL
    axes[i].spines["top"].set_visible(False)
    axes[i].spines["right"].set_visible(False)
    axes[i].legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        handlelength=1.5
    )

# Shared axis labels
fig.text(0.04, 0.5, "Normalized Absorbance", va="center", rotation="vertical", fontsize=font_size)
axes[-1].set_xlabel("Retention Volume (mL)")

#plt.subplots_adjust(hspace=0.37, left=0.08, right=0.95, bottom=0.07, top=0.97)
#plt.tight_layout(h_pad=1.0, rect=[0.08, 0.05, 1, 1])  # leave space for shared y label
#plt.tight_layout(h_pad=1.0, rect=[0.08, 0.05, 0.95, 0.97])
#plt.subplots_adjust(hspace=0.37)
plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.37)
plt.show()
