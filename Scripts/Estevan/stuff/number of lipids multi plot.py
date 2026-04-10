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
# Each entry corresponds to one trace on the plot.
# You can give each file its own protein mass and number of scaffold proteins.
datasets = [
    {
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\Control - Done\D1T0\20251105_EC_D1T0(-)1to80_DMPC_amberlite.1g_unidecfiles\20251105_EC_D1T0(-)1to80_DMPC_amberlite.1g_mass.txt",
        "label": "MSP1D1 DMPC",
        "protein_mass_per_scaffold_kDa": 22,  # mass of ONE scaffold protein for this dataset
        "n_scaffold_proteins": 2               # number of scaffold proteins per nanodisc for this dataset
    },
    {
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfGFP - done\D1T0\20251103_EC_D1sfGFP422uM_BCA_1to65_DMPC_native_unidecfiles\20251103_EC_D1sfGFP422uM_BCA_1to65_DMPC_native_mass.txt",
        "label": "MSP1D1-sfGFP DMPC",
        "protein_mass_per_scaffold_kDa": 48.8,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
{
        "file": r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Short - done\20251106_EC_D1sfCheryshort_1to80_DMPC_unidecfiles\20251106_EC_D1sfCheryshort_1to80_DMPC_mass.txt",
        "label": "MSP1D1-sfCherryS DMPC",
        "protein_mass_per_scaffold_kDa": 47.7,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
    # Add more entries here as needed:
    # {
    #     "file": r"path\to\another_mass.txt",
    #     "label": "Another protein / condition",
    #     "protein_mass_per_scaffold_kDa": 35.0,
    #     "n_scaffold_proteins": 2
    # },
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#A9A9A9",  # color for dataset 1
    "#009E60",  # color for dataset 2
    "#9B111E",  # extra colors if you add more datasets
    "#5C8D89",
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

handles = []
legend_labels = []

# Pre-calc lipid mass in kDa if needed
lipid_mass_kDa = lipid_mass_Da / 1000.0

# --- LOOP THROUGH DATASETS ---
for i, ds in enumerate(datasets):
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
        # Compute total protein mass per disc for THIS dataset (in Da)
        protein_total_mass_per_disc_Da = protein_mass_per_scaffold_kDa * 1000.0 * n_scaffold_proteins

        # number_of_lipids = (total_mass_Da - protein_total_mass_per_disc_Da) / lipid_mass_Da
        df["Num_Lipids"] = (df["Mass_Da"] - protein_total_mass_per_disc_Da) / lipid_mass_Da

        # Remove non-sensible points
        df = df[df["Num_Lipids"] > 0].copy()

        x_vals = df["Num_Lipids"]
        xlabel = "Number of lipids"
    else:
        x_vals = df["Mass_kDa"]
        xlabel = "Mass (kDa)"

    color = colors[i % len(colors)]

    # Plot
    line, = ax.plot(x_vals, df["Percent"], color=color, linewidth=1.5, label=label)
    handles.append(line)
    legend_labels.append(label)

# Final axes styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(x_range)
ax.set_ylim(0, 100)

ax.set_yticks([0, 50, 100])
ax.set_yticklabels(["0", "%", "100"])

ax.set_xlabel(xlabel, fontname=font_name, fontsize=font_size)
ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)

ax.legend(handles, legend_labels, loc="upper left", frameon=False, fontsize=font_size)

plt.tight_layout()
plt.show()