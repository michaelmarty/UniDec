import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

input_file = r"Z:\Group Share\Estevan\Lipid Stuffs\Molecule Transition Results.csv"

extract_sample = "HEK Extract"
discs_sample = "HEK_Discs_1"
# Optional:
# Enter an adduct such as "[M+H]" or "[M-H]" to analyze one mode only.
# Leave as None to include all precursor adducts.
adduct_to_use = None

# ------------------------------------------------------------
# Read the Skyline report
# ------------------------------------------------------------

df = pd.read_csv(input_file)

df["Area"] = pd.to_numeric(
    df["Area"],
    errors="coerce"
).fillna(0)

# ------------------------------------------------------------
# Keep precursor areas only
# ------------------------------------------------------------

df = df[
    df["Fragment Ion"].astype(str).str.lower().eq("precursor")
].copy()

# Keep only the two samples
df = df[
    df["Replicate Name"].isin(
        [extract_sample, discs_sample]
    )
].copy()

# Optionally keep only one precursor adduct
if adduct_to_use is not None:
    df = df[
        df["Precursor Adduct"] == adduct_to_use
    ].copy()

# ------------------------------------------------------------
# Calculate mTIC for each sample
# ------------------------------------------------------------

mtic = (
    df.groupby("Replicate Name")["Area"]
    .sum()
    .rename("mTIC")
)

df = df.merge(
    mtic,
    left_on="Replicate Name",
    right_index=True,
    how="left"
)

df["mTIC Normalized Area"] = (
    df["Area"] / df["mTIC"]
)

# ------------------------------------------------------------
# Sum normalized lipid areas by class
#
# "Molecule List Name" is the Skyline lipid-class column.
# ------------------------------------------------------------

class_data = (
    df.groupby(
        ["Molecule List Name", "Replicate Name"],
        as_index=False
    )["mTIC Normalized Area"]
    .sum()
)

# ------------------------------------------------------------
# Put the samples in separate columns
# ------------------------------------------------------------

comparison = class_data.pivot_table(
    index="Molecule List Name",
    columns="Replicate Name",
    values="mTIC Normalized Area",
    aggfunc="sum",
    fill_value=0
).reset_index()

for sample in [extract_sample, discs_sample]:
    if sample not in comparison.columns:
        comparison[sample] = 0

# Remove classes absent from either sample
comparison = comparison[
    (comparison[extract_sample] > 0) &
    (comparison[discs_sample] > 0)
].copy()

# ------------------------------------------------------------
# Calculate fold change
# ------------------------------------------------------------

comparison["Fold Change"] = (
    comparison[discs_sample] /
    comparison[extract_sample]
)

comparison["Log2 Fold Change"] = np.log2(
    comparison["Fold Change"]
)

comparison = comparison.sort_values(
    "Log2 Fold Change"
)

# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

print("\nmTIC values:")
print(mtic)

print("\nClass-level fold changes:")
print(
    comparison[
        [
            "Molecule List Name",
            extract_sample,
            discs_sample,
            "Fold Change",
            "Log2 Fold Change"
        ]
    ].to_string(index=False)
)

# ------------------------------------------------------------
# Plot class-level fold changes
# ------------------------------------------------------------


plt.figure(figsize=(7, 4))

plt.bar(
    comparison["Molecule List Name"],
    comparison["Log2 Fold Change"],
    width=0.6
)

plt.axhline(0, linewidth=1)

plt.xlabel("Lipid Class")
plt.ylabel("Log2 Fold Change")
plt.title("mTIC-Normalized Lipid Class Fold Changes")

plt.xticks(rotation=45, ha="right", fontsize=8)
plt.yticks(fontsize=8)

plt.tight_layout()
plt.show()