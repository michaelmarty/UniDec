import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10

# If x_as_lipids is True, the x-axis will be converted to number of lipids.
x_as_lipids = True

# If using lipid-count x-axis, provide these:
# Protein mass per scaffold protein in kDa (enter the mass of ONE scaffold protein)
protein_mass_per_scaffold_kDa = 30.2   # example: 25 kDa (change to your protein)
# Number of scaffold proteins per nanodisc (often 2 for MSP-based discs)
n_scaffold_proteins = 2
# Lipid mass in Daltons (DMPC ≈ 677.93 Da). Enter as Da.
lipid_mass_Da = 678

# If you still want to set x-axis limits (in number of lipids when x_as_lipids=True,
# otherwise in kDa), set x_range appropriately:
# Example for lipid counts:
x_range = (250, 450)  # number of lipids limits (if x_as_lipids True)
# Example for mass (kDa): x_range = (140, 190)

# === FILES TO PLOT ===
files = [
    r"Z:\Group Share\Estevan\Q Exactive Data\20251205\20251205_EC_E3t3_1to150DMPC_unidecfiles\20251205_EC_E3t3_1to150DMPC_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251205\20251205_EC_E3t3_1to150DMPC_cdms_2025-12-05-11-07-52_unidecfiles\20251205_EC_E3t3_1to150DMPC_cdms_2025-12-05-11-07-52_mass.txt"
]

# === LEGEND LABELS (one per file) ===
labels = [
    "MSP1E3D1 DMPC Native-MS",
    "MSP1E3D1 DMPC CD-MS"
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#CACCC8",  # color for file 1
    "#85949F",  # color for file 2
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

# Keep handles for the legend
handles = []
legend_labels = []

# Pre-calc constants and convert protein mass to Da
protein_total_mass_per_disc_Da = protein_mass_per_scaffold_kDa * 1000.0 * n_scaffold_proteins
lipid_mass_kDa = lipid_mass_Da / 1000.0

# --- LOOP THROUGH FILES ---
for i, (file, label) in enumerate(zip(files, labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert mass to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    if x_as_lipids:
        # Compute number of lipids for each mass point:
        # number_of_lipids = (total_mass_Da - protein_total_mass_per_disc_Da) / lipid_mass_Da
        df["Num_Lipids"] = (df["Mass_Da"] - protein_total_mass_per_disc_Da) / lipid_mass_Da

        # Remove points where the computed number_of_lipids is non-sensible (<=0)
        df = df[df["Num_Lipids"] > 0].copy()

        x_vals = df["Num_Lipids"]
        xlabel = "Number of lipids"
    else:
        x_vals = df["Mass_kDa"]
        xlabel = "Mass (kDa)"

    color = colors[i % len(colors)]

    # Plot the main curve on the single axes and capture handle for legend
    line, = ax.plot(x_vals, df["Percent"], color=color, linewidth=1.5, label=label)
    handles.append(line)
    legend_labels.append(label)

# Final axes styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)
ax.set_ylim(0, 100)           # ensure 0-100 percent range

# Set desired ticks and labels: ticks at 0, 50, 100 and labels "0", "%", "100"
ax.set_yticks([0, 50, 100])
ax.set_yticklabels(["0", "%", "100"])

ax.set_xlabel(xlabel, fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

# Single legend (shows a line with matching color for each trace)
ax.legend(handles, legend_labels, loc="upper left", frameon=False, fontsize=font_size)

plt.tight_layout()
plt.show()