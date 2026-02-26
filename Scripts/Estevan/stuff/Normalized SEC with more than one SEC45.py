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

# -----------------------
# FILES TO PLOT (TRIPLICATES PER SAMPLE)
# -----------------------
files = [
    [
        r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to80DMPC.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to72_DMPC.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1t65_DMPC.xlsx",
    ],
    [
        r"D:\UT Spring 2026\general lab data\XLSX akta files\E3sfGFP_DMPC1.xlsx",
        r"D:\UT Spring 2026\general lab data\XLSX akta files\E3sfGFP_DMPC_3.xlsx",
        r"D:\UT Spring 2026\general lab data\XLSX akta files\E3sfGFP_DMPC_4.xlsx",
    ],
    # add more samples here
]

# -----------------------
# LEGEND LABELS (ONE PER SUBPLOT)
# -----------------------
sample_labels = [
    "MSP1D1-sfGFP",
    "MSP1E3D1-sfGFP",
]

# -----------------------
# COLORS (ONE PER SAMPLE)
# -----------------------
colors = [
    "#333333",
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

for i, (replicate_files, label) in enumerate(zip(files, sample_labels)):

    for file in replicate_files:
        df = read_sec_file(file)
        df["mAU"] = df["mAU"] / df["mAU"].max()

        # Line
        axes[i].plot(
            df["ml"],
            df["mAU"],
            color=colors[i],
            linewidth=1.4,
            alpha=0.8
        )

        # Fill
        axes[i].fill_between(
            df["ml"],
            df["mAU"],
            0,
            color=colors[i],
            alpha=0.15,
            linewidth=0
        )

    # -----------------------
    # Axis formatting
    # -----------------------
    axes[i].set_ylim(0, 1.05)
    axes[i].margins(y=0)
    axes[i].set_xlim(10, 20)

    axes[i].spines["top"].set_visible(False)
    axes[i].spines["right"].set_visible(False)

    axes[i].set_yticks([0, 1])
    axes[i].tick_params(axis="y", length=3, width=0.8)
    axes[i].tick_params(axis="x", length=3, width=0.8)

    axes[i].legend(
        [label],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        handlelength=1.5
    )

# -----------------------
# Shared axis labels
# -----------------------
fig.text(
    0.04, 0.5,
    "Normalized Absorbance",
    va="center",
    rotation="vertical",
    fontsize=font_size
)

axes[-1].set_xlabel("Retention Volume (mL)")

# -----------------------
# Layout (UNCHANGED)
# -----------------------
plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)

plt.show()
