import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
x_range = (150, 250)  # x-axis limits (in kDa)

# === FILES TO PLOT ===
files = [
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Short - done\20251106_EC_MSP1D1sfCherryS_CDMS_rep1_unidecfiles\20251106_EC_MSP1D1sfCherryS_CDMS_rep1_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Short - done\20251106_EC_MSP1D1sfCherryS_NATIVE_rep1_unidecfiles\20251106_EC_MSP1D1sfCherryS_NATIVE_rep1_mass.txt"
]

# === LEGEND LABELS (one per file) ===
labels = [
    "MSP1D1sfCheryS CDMS",
    "MSP1D1sfCherryS Native MS"
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#A16084",  # color for file 1
    "#CC79A7",  # color for file 2
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# Keep handles for the legend
handles = []
legend_labels = []

# --- LOOP THROUGH FILES ---
for i, (file, label) in enumerate(zip(files, labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    color = colors[i % len(colors)]

    # Plot the main curve on the single axes and capture handle for legend
    line, = ax.plot(df["Mass_kDa"], df["Percent"], color=color, linewidth=1.5, label=label)
    handles.append(line)
    legend_labels.append(label)

# Final axes styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)          # kDa limits
ax.set_ylim(0, 100)           # ensure 0-100 percent range

# Set desired ticks and labels: ticks at 0, 50, 100 and labels "0", "%", "100"
ax.set_yticks([0, 50, 100])
ax.set_yticklabels(["0", "%", "100"])

ax.set_xlabel("Mass (kDa)", fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

# Single legend (shows a line with matching color for each trace)
ax.legend(handles, legend_labels, loc="upper left", frameon=False, fontsize=font_size)

plt.tight_layout()
plt.show()