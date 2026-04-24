import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10

# If x_as_lipids is True, the x-axis will be converted to number of lipids.
x_as_lipids = True

# === PER-FILE PROTEIN PARAMETERS ===
# One entry per file, in the same order as `files`.
protein_mass_per_scaffold_kDa_list = [
    22,   # for file 1
    47.7,   # for file 2
    48.8,   # for file 3 (example; add/remove as needed)
]

n_scaffold_proteins_list = [
    2,  # for file 1
    2,  # for file 2
    2,  # for file 3
]

# Lipid mass in Daltons (same for all files; change if needed)
lipid_mass_Da = 678

# x-axis limits (lipid counts if x_as_lipids=True, otherwise in kDa)
x_range = (100, 240)

# === FILES TO PLOT ===
files = [
     r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\Control - Done\D1T0\20251105_EC_MSP1D1_DMPC_CDMS_rep2_unidecfiles\20251105_EC_MSP1D1_DMPC_CDMS_rep2_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfCherry - done\Short - done\20251106_EC_MSP1D1sfCherryS_CDMS_rep1_unidecfiles\20251106_EC_MSP1D1sfCherryS_CDMS_rep1_mass.txt",
    r"Z:\Group Share\Estevan\Fluoro_ND_paper\Replicates\FluroND\sfGFP - done\D1T0\20251203EC_MSP1D1sfGFP_CDMS_rep3_unidecfiles\20251203EC_MSP1D1sfGFP_CDMS_rep3_mass.txt",  # example third file
]

# === LEGEND LABELS (one per file) ===
labels = [
    "MSP1D1",
    "MSP1D1-sfCherryS",
    "MSP1D1-sfGFP",
]

# Sanity check: lengths must match
assert len(files) == len(labels) == len(protein_mass_per_scaffold_kDa_list) == len(n_scaffold_proteins_list), \
    "files, labels, and protein parameter lists must have the same length."

n_files = len(files)

# === FIGURE SETUP: N stacked subplots sharing x-axis ===
fig, axes = plt.subplots(n_files, 1, figsize=(6, 2.2 * n_files), dpi=300, sharex=True)

# If only one file, axes is a single Axes object; make it a list for consistency
if n_files == 1:
    axes = [axes]

# === CUSTOM COLORS (repeat if fewer than files) ===
colors = [
    "#000000",
    "#CF7CA9",
    "#00A9A5",
]
# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

xlabel = ""

# --- LOOP THROUGH FILES AND PLOT ON SEPARATE AXES ---
for i, (file, label) in enumerate(zip(files, labels)):
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert mass to kDa and normalize
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    # Per-file protein parameters
    protein_mass_per_scaffold_kDa = protein_mass_per_scaffold_kDa_list[i]
    n_scaffold_proteins = n_scaffold_proteins_list[i]

    # Total protein mass per disc in Da for this file
    protein_total_mass_per_disc_Da = protein_mass_per_scaffold_kDa * 1000.0 * n_scaffold_proteins

    if x_as_lipids:
        # Compute number of lipids for each mass point
        df["Num_Lipids"] = (df["Mass_Da"] - protein_total_mass_per_disc_Da) / lipid_mass_Da
        df = df[df["Num_Lipids"] > 0].copy()

        x_vals = df["Num_Lipids"]
        xlabel = "Number of lipids"
    else:
        x_vals = df["Mass_kDa"]
        xlabel = "Mass (kDa)"

    color = colors[i % len(colors)]
    ax = axes[i]

    ax.plot(x_vals, df["Percent"], color=color, linewidth=2, label=label)

    # Styling for each subplot
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(["0", "%", "100"])
    ax.set_ylabel("Relative Intensity (%)", fontname=font_name, fontsize=font_size)
    ax.legend(loc="upper left", frameon=False, fontsize=font_size)

# Shared x-axis limits on the bottom axis
axes[-1].set_xlim(x_range)
axes[-1].set_xlabel(xlabel, fontname=font_name, fontsize=font_size)

# Hide x tick labels on all but the bottom axis
for ax in axes[:-1]:
    plt.setp(ax.get_xticklabels(), visible=False)

plt.tight_layout()
plt.show()