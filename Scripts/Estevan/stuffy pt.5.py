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
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1T0minus_1to80_DMPC 001.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to80DMPC.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 003sfCherryshort_1to80_DMPC_rep2 002.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 003_D1blueberry_1to80_DMPC_BCAchecked 001.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1halo1toDMPC 001.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_BIORAD_secstdtry1 001.xlsx",
]

# === LEGEND LABELS ===
sample_labels = [
    "MSP1D1",
    "MSP1D1-sfGFP",
    "MSP1D1-sfCherryshort",
    "MSP1D1-mBlueberry",
    "MSP1D1-HaloTag",
    #"Standards"
]

# === CUSTOM HEX COLORS ===
colors = [
    "#333333",  # charcoal
    "#00A9A5",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # blue
    "#999999",  # grey
    #"#4D4D4D"   # yellow-green
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
    axes[i].set_xlim(0, 25)  # cutoff at 25 mL
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
