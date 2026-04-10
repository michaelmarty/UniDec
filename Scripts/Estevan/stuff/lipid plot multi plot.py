import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10

# If x_as_lipids is True, the x-axis will be converted to number of lipids.
x_as_lipids = True

# Lipid mass in Daltons (DMPC ≈ 677.93 Da). Enter as Da.
lipid_mass_Da = 678

# If you still want to set x-axis limits (in number of lipids when x_as_lipids=True,
# otherwise in kDa), set x_range appropriately:
# Example for lipid counts:
x_range = (100, 200)  # number of lipids limits (if x_as_lipids True)
# Example for mass (kDa): x_range = (140, 190)

# === FILES & PROTEIN PARAMETERS ===
datasets = [
    {
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\Control - Done\D1T0\20251105_EC_D1T0(-)1to80_DMPC_amberlite.1g_unidecfiles\20251105_EC_D1T0(-)1to80_DMPC_amberlite.1g_mass.txt",
        "label": "MSP1D1 DMPC",
        "protein_mass_per_scaffold_kDa": 22,  # mass of ONE scaffold protein for this dataset
        "n_scaffold_proteins": 2  # number of scaffold proteins per nanodisc for this dataset
    },
    {
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfGFP - done\D1T0\20251103_EC_D1sfGFP422uM_BCA_1to65_DMPC_native_unidecfiles\20251103_EC_D1sfGFP422uM_BCA_1to65_DMPC_native_mass.txt",
        "label": "MSP1D1-sfGFP DMPC",
        "protein_mass_per_scaffold_kDa": 48.8,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
    {
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Short - done\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_unidecfiles\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_mass.txt",
        "label": "MSP1D1-sfCherryS DMPC",
        "protein_mass_per_scaffold_kDa": 47.7,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
    # Add more datasets as needed...
]

# === SUBPLOT LAYOUT ===
# For N datasets, you can choose:
# - nrows = N, ncols = 1  (stacked vertically)
# - or something like nrows = 2, ncols = 2 if you have 4 datasets, etc.
n_datasets = len(datasets)
nrows = n_datasets
ncols = 1

# === FIGURE SETUP ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 3 * nrows), dpi=300, sharex=True)

# If there's only one dataset, axes is not an array; make it iterable
if n_datasets == 1:
    axes = [axes]

# === CUSTOM COLORS ===
colors = [
    "#A9A9A9",  # color for dataset 1
    "#009E60",  # color for dataset 2
    "#9B111E",
    "#5C8D89",
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

lipid_mass_kDa = lipid_mass_Da / 1000.0

# To ensure consistent x-axis label text
xlabel = "Number of lipids" if x_as_lipids else "Mass (kDa)"

# --- LOOP THROUGH DATASETS AND AXES ---
for i, (ds, ax) in enumerate(zip(datasets, axes)):
    file = ds["file"]
    label = ds["label"]
    protein_mass_per_scaffold_kDa = ds["protein_mass_per_scaffold_kDa"]
    n_scaffold_proteins = ds["n_scaffold_proteins"]

    # Read file
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert mass to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    if x_as_lipids:
        protein_total_mass_per_disc_Da = protein_mass_per_scaffold_kDa * 1000.0 * n_scaffold_proteins
        df["Num_Lipids"] = (df["Mass_Da"] - protein_total_mass_per_disc_Da) / lipid_mass_Da
        df = df[df["Num_Lipids"] > 0].copy()
        x_vals = df["Num_Lipids"]
    else:
        x_vals = df["Mass_kDa"]

    color = colors[i % len(colors)]

    # Plot on this axis
    ax.plot(x_vals, df["Percent"], color=color, linewidth=1.5, label=label)

    # Styling for this axis
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(["0", "%", "100"])
    ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

    # Legend for this subplot
    ax.legend(loc="upper left", frameon=False, fontsize=font_size)

# Apply x-axis limits and label on the bottom subplot
axes[-1].set_xlim(x_range)
axes[-1].set_xlabel(xlabel, fontname=font_name, fontsize=font_size)

plt.tight_layout()
plt.show()