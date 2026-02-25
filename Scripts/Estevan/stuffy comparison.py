import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === FILES TO PLOT (GROUP REPLICATES TOGETHER) ===
files =[
    [r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt", r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_Imidazole_2025-09-02-01-41-53_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_Imidazole_2025-09-02-01-41-53_mass.txt"],
    [r"Z:\Group Share\Estevan\Q Exactive Data\20251020\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_20251020154408_2025-10-20-04-12-28_unidecfiles\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_20251020154408_2025-10-20-04-12-28_mass.txt", r"Z:\Group Share\Estevan\Q Exactive Data\20251020\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_maybenumba3_2025-10-20-04-46-45_unidecfiles\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_maybenumba3_2025-10-20-04-46-45_mass.txt"],
    [r"Z:\Group Share\Estevan\Q Exactive Data\20251022\20251022_D1sfGFPminus_1to150_DMPC_CDMS_2025-10-22-02-34-17_unidecfiles\20251022_D1sfGFPminus_1to150_DMPC_CDMS_2025-10-22-02-34-17_mass.txt", r"Z:\Group Share\Estevan\Q Exactive Data\20251020\20251020_EC_D1sfGFPminustag1to150_CDMS_1to3dilution_2025-10-20-04-11-35_unidecfiles\20251020_EC_D1sfGFPminustag1to150_CDMS_1to3dilution_2025-10-20-04-11-35_mass.txt"]
]
# === LEGEND LABELS (ONE PER CONDITION) ===
labels = ["MSP1D1",
    "MSP1D1-sfCherry",
"MSP1D1-sfGFP",
    # "Sample A",
    # "Sample B",
    # "Sample C"
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
    "#440154FF",  # deep purple
    "#414487FF",  # dark blue
    "#2A788EFF"   # blue-green
]

# === Y-OFFSET SETTINGS ===
y_offset_step = 120
y_labels = ["0", "%", "100"]

# === LOOP THROUGH CONDITIONS ===
for i, (replicate_group, label) in enumerate(zip(files, labels)):
    ax = axes[i]
    offset = i * y_offset_step
    color = colors[i % len(colors)]

    # === PLOT EACH REPLICATE ===
    for j, file in enumerate(replicate_group):
        df = pd.read_csv(file, sep=r"\s+", header=None)
        df.columns = ["Mass_Da", "Intensity"]

        # Convert to kDa and normalize to 100%
        df["Mass_kDa"] = df["Mass_Da"] / 1000
        df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

        alpha_val = 0.6 if j < len(replicate_group) - 1 else 1.0  # last replicate = fully opaque
        ax.plot(
            df["Mass_kDa"],
            df["Percent"] + offset,
            color=color,
            linewidth=1.5,
            alpha=alpha_val
        )

    # === STYLE ===
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(50, 350)
    ax.set_ylim(offset, offset + 100)
    ax.set_yticks([offset, offset + 50, offset + 100])
    ax.set_yticklabels(y_labels)

    # === LEGEND (ONE PER CONDITION) ===
    ax.legend([label], loc="center left", bbox_to_anchor=(1.15, 0.5), frameon=False)

    # Only show x-axis labels on bottom plot
    if i < len(files) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Mass (kDa)")

# === COMMON Y-AXIS LABEL ===
fig.text(0.03, 0.5, "Relative Intensity (%)", va="center", rotation="vertical")

# === TRANSPARENT BACKGROUND ===
fig.patch.set_alpha(0.0)
for ax in axes:
    ax.set_facecolor("none")

# === FIX CLIPPING & SPACING ===
plt.subplots_adjust(left=0.2, right=0.88, bottom=0.15, top=0.95, hspace=0.45)
plt.tight_layout(rect=[0.05, 0, 0.9, 1])

plt.show()
