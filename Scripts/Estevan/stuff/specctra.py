import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === FILES TO PLOT ===
files = [
      r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251020\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_20251020154408_2025-10-20-04-12-28_unidecfiles\20251020_EC_D1sfCherryshortminustag1to150_CDMS_1to3dilution_20251020154408_2025-10-20-04-12-28_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251020\20251020_EC_D1sfGFPminustag1to150_CDMS_1to3dilution_2025-10-20-04-11-35_unidecfiles\20251020_EC_D1sfGFPminustag1to150_CDMS_1to3dilution_2025-10-20-04-11-35_mass.txt",
]

# === FIGURE SETUP ===
fig, axes = plt.subplots(len(files), 1, sharex=True, figsize=(5, 4), dpi=300)
if len(files) == 1:
    axes = [axes]

# === Y-OFFSET SETTINGS ===
y_offset_step = 120  # vertical spacing between traces
y_labels = ["0", "%", "100"]

# === LOOP THROUGH FILES ===
for i, file in enumerate(files):
    print(f"📄 Reading file: {file}")
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize to 100%
    df["Mass_kDa"] = df["Mass_Da"] / 1000
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    ax = axes[i]
    offset = i * y_offset_step
    ax.plot(df["Mass_kDa"], df["Percent"] + offset, color="black", linewidth=1.2)

    # Axis styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 300)
    ax.set_ylim(offset, offset + 100)
    ax.set_yticks([offset, offset + 50, offset + 100])
    ax.set_yticklabels(y_labels)

    # Only show x-axis tick labels on the bottom plot
    if i < len(files) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Mass (kDa)")

# Common y-axis label
fig.text(0.01, 0.5, "Relative Intensity (%)", va="center", rotation="vertical")

# Make background transparent
fig.patch.set_alpha(0.0)
for ax in axes:
    ax.set_facecolor("none")

# Adjust layout: more padding for left label and spacing between plots
plt.subplots_adjust(left=0.15, hspace=0.4, bottom=0.12)

plt.show()
