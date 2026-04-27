import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
x_range = (58, 60)  # x-axis limits (in kDa)

# === FILES TO PLOT ===
files = [
    r"C:\Python\UniDecPastedSpectra\20260323_MSP1E3D1_mCherry_plustag_unidecfiles\20260323_MSP1E3D1_mCherry_plustag_mass.txt",
   # r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\Mass txt files\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_mass.txt"
]

# === MAIN LEGEND LABELS (one per file) ===
labels = [
    "MSP1E3D1-mCherry(+)",
    #"MSP1D1-mCherry"
]

# === SECOND LEGEND LABELS (one per file) ===
second_labels = [
    "MSP1D1-sfGFp DMPC + AM2",
   # "151 ± 9 DMPC per Nanodisc",
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(9, 5), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#5E2F98",  # color for file 1
    #"#0072B2",  # color for file 2
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# --- LOOP THROUGH FILES ---
for i, (file, label, second_label) in enumerate(zip(files, labels, second_labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    color = colors[i % len(colors)]

    # Plot the main curve on the single axes (no legend)
    ax.plot(df["Mass_kDa"], df["Percent"], color=color, linewidth=1.5)

# Final axes styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)  # kDa limits
ax.set_ylim(bottom=0)

ax.set_xlabel("Mass (kDa)", fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

plt.tight_layout()
plt.show()