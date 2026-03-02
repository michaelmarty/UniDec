import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
x_range = (0, 500)  # x-axis limits

# spacing_factor and show_per_trace_ticks are no longer needed for Option 1
# spacing_factor = 0.5
# show_per_trace_ticks = True

# === FILES TO PLOT ===
files = [
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfCherry - done\Long - done\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_unidecfiles\20251208_EC_D1cherryLong_1to90_DMPC_CDMS_2025-12-08-11-04-12_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mCherry - done\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_unidecfiles\20251202_EC_D1mCherryREDOcdms_2025-12-02-11-15-03_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\mBlueberry - done\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_unidecfiles\20251107_EC_D1blueberry_1to80_DMPC_CDMS_2025-11-07-12-07-14_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\Halo-Tag - done\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_unidecfiles\20251204_EC_D1Halo_1to80_DMPC_rep3cdms_2025-12-04-11-13-01_mass.txt",
r"Z:\Group Share\Estevan\Aiden\Replicates\FluroND\sfGFP - done\E3T0\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_unidecfiles\20251209_EC_E3sfGFP_1to160_DMPC_run_twoCDMS_2025-12-09-01-33-49_mass.txt"
]

# === MAIN LEGEND LABELS (one per file) ===
labels = [
    "MSP1D1-sfCherryL",
    "MSP1D1-mCherry",
    "MSP1D1-mBlueberry",
    "MSP1D1-HaloTag",
    "MSP1E3D1-sfGFP"
]

# === SECOND LEGEND LABELS (one per file) ===
second_labels = [
    "185 ± 16 DMPC per Nanodisc",
    "154 ± 9 DMPC per Nanodisc",
    "156 ± 7 DMPC per Nanodisc",
    "146 ± 14 DMPC per Nanodisc",
    "360 ± 20 DMPC per Nanodisc"
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(9, 5), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#CC79A7",
    "#CC79A7",
    "#2C7FB8",
    "#7A7A7A",
    "#00A9A5"
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# Read all data, normalize intensity to percent
data_list = []
for file in files:
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100.0
    data_list.append(df)

# Plot traces with NO vertical offset (shared 0–100% y-axis)
handles = []
for i, (df, label, second_label) in enumerate(zip(data_list, labels, second_labels)):
    color = colors[i % len(colors)]
    h, = ax.plot(df["Mass_kDa"], df["Percent"], color=color, linewidth=1.5)
    handles.append(h)

    # annotate DMPC label near the peak
    peak_idx = df["Percent"].idxmax()
    peak_mass = df["Mass_kDa"].iloc[peak_idx]
    peak_value = df["Percent"].iloc[peak_idx]
    x_min, x_max = x_range
    x_shift = (x_max - x_min) * 0.02
    ax.text(peak_mass + x_shift, peak_value, second_label,
            fontsize=font_size * 0.85, color=color, va='center', ha='left', alpha=0.95)

# Single legend for traces
ax.legend(handles, labels, loc='upper right', frameon=False, fontsize=font_size)

# Axis styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)

# Y limits and ticks: common 0–100% scale
ax.set_ylim(-5, 105)  # small padding
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)

# Remove custom per-trace ticks (no longer needed)
# (You can delete the whole show_per_trace_ticks block from your original code.)

# Axis labels
ax.set_xlabel("Mass (kDa)", fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)
plt.show()