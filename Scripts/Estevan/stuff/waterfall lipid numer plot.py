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

# X-axis limits:
# If x_as_lipids=True, this is in number of lipids; otherwise in kDa.
x_range = (100, 200)

# Vertical offset between traces (in "percent intensity" units)
offset = -2  # try 20–40 depending on how many traces you have

# === FILES & PROTEIN PARAMETERS ===
datasets = [
    {
        "file": r"C:\Users\ec45294\The University of Texas at Austin\MartyLab - General\Presentations\Undergrads\Figures\Data\20251105_EC_MSP1D1_DMPC_NATIVErep1_unidecfiles\20251105_EC_MSP1D1_DMPC_NATIVErep1_mass.txt",
        "label": "MSP1D1 DMPC",
        "protein_mass_per_scaffold_kDa": 22,  # mass of ONE scaffold protein for this dataset
        "n_scaffold_proteins": 2  # number of scaffold proteins per nanodisc for this dataset
    },
    {
        "file": r"C:\Users\ec45294\The University of Texas at Austin\MartyLab - General\Presentations\Undergrads\Figures\Data\20251106_EC_MSP1D1sfCherryS_NATIVE_rep1_unidecfiles\20251106_EC_MSP1D1sfCherryS_NATIVE_rep1_mass.txt",
        "label": "MSP1D1-sfCherryS DMPC",
        "protein_mass_per_scaffold_kDa": 47.7,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
    {
        "file": r"C:\Users\ec45294\The University of Texas at Austin\MartyLab - General\Presentations\Undergrads\Figures\Data\20251103_EC_MSP1D1sfGFP_NATIVE_rep1_unidecfiles\20251103_EC_MSP1D1sfGFP_NATIVE_rep1_mass.txt",
        "label": "MSP1D1-sfGFP DMPC",
        "protein_mass_per_scaffold_kDa": 48.8,  # example: different protein or construct
        "n_scaffold_proteins": 2
    },
    # Add more datasets as needed...
]

# === FIGURE SETUP ===
fig, ax = plt.subplots(1, 1, figsize=(7, 8), dpi=300)

# === CUSTOM COLORS ===
colors = [
    "#6e6e77",  # color for dataset 1
    "#b34297",  # color for dataset 2
    "#159d28",
    "#5C8D89",
]

# Apply font settings globally
plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

lipid_mass_kDa = lipid_mass_Da / 1000.0
xlabel = "Number of lipids" if x_as_lipids else "Mass (kDa)"

handles = []
legend_labels = []

# --- LOOP THROUGH DATASETS ---
for i, ds in enumerate(datasets):
    file = ds["file"]
    label = ds["label"]
    protein_mass_per_scaffold_kDa = ds["protein_mass_per_scaffold_kDa"]
    n_scaffold_proteins = ds["n_scaffold_proteins"]

    # Read file
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df.columns = ["Mass_Da", "Intensity"]

    # Convert mass to kDa and normalize to 0–100 %
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

    if x_as_lipids:
        protein_total_mass_per_disc_Da = (
            protein_mass_per_scaffold_kDa * 1000.0 * n_scaffold_proteins
        )
        df["Num_Lipids"] = (df["Mass_Da"] - protein_total_mass_per_disc_Da) / lipid_mass_Da
        df = df[df["Num_Lipids"] > 0].copy()
        x_vals = df["Num_Lipids"]
    else:
        x_vals = df["Mass_kDa"]

    # Apply vertical offset for waterfall effect
    y_vals = df["Percent"] + i * offset

    color = colors[i % len(colors)]

    line, = ax.plot(x_vals, y_vals, color=color, linewidth=1.5, label=label)
    handles.append(line)
    legend_labels.append(label)

# Styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(x_range)

# Y‑axis limits: from slightly below the first trace to above the last
n_traces = len(datasets)
ax.set_ylim(-5, 100 + (n_traces - 1) * offset + 5)

# Optional: show y‑ticks at each trace baseline (0, offset, 2*offset,...)
# or just leave default. Here we keep it simple:
ax.set_ylabel("Relative Intensity (offset waterfall)", fontname=font_name, fontsize=font_size)

ax.set_xlabel(xlabel, fontname=font_name, fontsize=font_size)

ax.legend(handles, legend_labels, loc="upper left", frameon=False, fontsize=font_size)

plt.tight_layout()
plt.show()