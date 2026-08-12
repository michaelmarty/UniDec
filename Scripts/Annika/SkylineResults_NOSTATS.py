import pandas as pd

# -----------------------------
# Input and output files
# -----------------------------
input_file = r"Z:\Group Share\Annika\Stellar\TAMU\Synthetic - Basic\Area Results_Simplified.csv"
output_file = r"Z:\Group Share\Annika\Stellar\TAMU\Synthetic - Basic\Area Results_Normalized_Percents.csv"

# -----------------------------
# Read CSV
# -----------------------------
df = pd.read_csv(input_file)

# -----------------------------
# Clean Area column
# Turns #N/A, N/A, blanks, text errors, etc. into 0
# -----------------------------
df["Area"] = pd.to_numeric(df["Area"], errors="coerce").fillna(0)

# -----------------------------
# Get lipid class from Molecule
# Examples:
# IS CL 14:0-14:0-14:0-14:0 -> CL
# CL 16:0-16:0-16:0-16:0    -> CL
# IS PA 15:0-18:1(+[2]H7)    -> PA
# PA 16:0-18:1               -> PA
# -----------------------------
def get_class(molecule):
    parts = str(molecule).split()

    if parts[0] == "IS":
        return parts[1]
    else:
        return parts[0]

df["Class"] = df["Molecule"].apply(get_class)

# -----------------------------
# Mark internal standards
# -----------------------------
df["Is_IS"] = df["Molecule"].astype(str).str.startswith("IS ")

# -----------------------------
# Combine duplicate rows
# Same Molecule + Replicate gets summed
# -----------------------------
combined = (
    df.groupby(["Molecule", "Replicate", "Class", "Is_IS"], as_index=False)["Area"]
    .sum()
)

# -----------------------------
# Make IS lookup table
# One IS per Class and Replicate
# -----------------------------
is_table = combined[combined["Is_IS"]].copy()

is_table = is_table.rename(columns={
    "Molecule": "IS_Molecule",
    "Area": "IS_Area"
})

is_table = is_table[["Class", "Replicate", "IS_Molecule", "IS_Area"]]

# -----------------------------
# Attach matching IS area to each lipid
# -----------------------------
normalized = combined.merge(
    is_table,
    on=["Class", "Replicate"],
    how="left"
)

# -----------------------------
# IS normalize
# -----------------------------
normalized["IS_Normalized_Area"] = (
    normalized["Area"] / normalized["IS_Area"]
)

# If no IS or IS area is 0, set normalized value to 0
normalized.loc[
    normalized["IS_Area"].isna() | (normalized["IS_Area"] == 0),
    "IS_Normalized_Area"
] = 0

# -----------------------------
# Remove internal standards from final output
# -----------------------------
normalized = normalized[~normalized["Is_IS"]].copy()

# -----------------------------
# Save all IS-normalized lipids
# -----------------------------
normalized = normalized[
    [
        "Molecule",
        "Class",
        "Replicate",
        "Area",
        "IS_Molecule",
        "IS_Area",
        "IS_Normalized_Area"
    ]
]

normalized.to_csv(output_file, index=False)

# -----------------------------
# Optional: calculate percent of selected lipids only
# Example:
# PC / (PC + CL + PA) * 100
# -----------------------------
lipids_to_include = [
    "PA 16:0-18:1",
    "CL 16:0-16:0-16:0-16:0",
    "PC 16:0-18:1"
]

selected = normalized[
    normalized["Molecule"].isin(lipids_to_include)
].copy()

selected_total = (
    selected.groupby("Replicate", as_index=False)["IS_Normalized_Area"]
    .sum()
    .rename(columns={"IS_Normalized_Area": "Selected_Total_IS_Normalized"})
)

selected = selected.merge(
    selected_total,
    on="Replicate",
    how="left"
)

selected["Percent_of_Selected_Total"] = (
    selected["IS_Normalized_Area"] /
    selected["Selected_Total_IS_Normalized"] *
    100
)

selected.loc[
    selected["Selected_Total_IS_Normalized"].isna() |
    (selected["Selected_Total_IS_Normalized"] == 0),
    "Percent_of_Selected_Total"
] = 0

selected = selected[
    [
        "Molecule",
        "Class",
        "Replicate",
        "Area",
        "IS_Molecule",
        "IS_Area",
        "IS_Normalized_Area",
        "Selected_Total_IS_Normalized",
        "Percent_of_Selected_Total"
    ]
]

selected.to_csv(output_file, index=False)

print("Done!")
print(f"Saved IS-normalized output as: {output_file}")
print(f"Saved selected lipid percent output as: {output_file}")