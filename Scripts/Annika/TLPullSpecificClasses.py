import pandas as pd
import numpy as np

# ===== USER INPUT =====
input_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Combined_TransitionList_Cleanest_2.csv"
output_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\PCs-SMs_TransList.csv"

MOLECULE_COLUMN = "Molecule List Name"
ADDUCT_COLUMN = "Precursor Adduct"
MZ_COLUMN = "Precursor Mz"
PRECURSOR_CHARGE_COLUMN = "Precursor Charge"
PRODUCT_CHARGE_COLUMN = "Product Charge"

classes_to_keep = ["PC", "LPC", "SM"]

PROTON = 1.007276
ACETATE = 59.013304
MZ_TOL = 0.0005
# ======================

df = pd.read_csv(input_csv)

# Filter classes
filtered_df = df[df[MOLECULE_COLUMN].isin(classes_to_keep)].copy()

# Split positives and existing acetates
positive_df = filtered_df[
    filtered_df[ADDUCT_COLUMN] == "[M+H]+"
].copy()

existing_acetate = filtered_df[
    filtered_df[ADDUCT_COLUMN] == "[M+CH3COO]-"
].copy()

new_rows = []

# Group by precursor (molecule + precursor m/z)
group_cols = ["Molecule Name", MZ_COLUMN]

for (mol_name, precursor_mz), group in positive_df.groupby(group_cols):

    # Calculate new acetate m/z
    neutral_mass = precursor_mz - PROTON
    new_mz = round(neutral_mass + ACETATE, 5)

    # Check if acetate already exists for this molecule
    acetate_match = existing_acetate[
        (existing_acetate["Molecule Name"] == mol_name) &
        (np.abs(existing_acetate[MZ_COLUMN] - new_mz) < MZ_TOL)
    ]

    if len(acetate_match) == 0:

        # Duplicate entire precursor group
        new_group = group.copy()

        new_group[MZ_COLUMN] = new_mz
        new_group[ADDUCT_COLUMN] = "[M+CH3COO]-"
        new_group[PRECURSOR_CHARGE_COLUMN] = -1
        new_group[PRODUCT_CHARGE_COLUMN] = -1

        new_rows.append(new_group)

# Append new acetate transitions
if len(new_rows) > 0:
    acetate_df = pd.concat(new_rows, ignore_index=True)
    final_df = pd.concat([filtered_df, acetate_df], ignore_index=True)
else:
    final_df = filtered_df

# Optional rounding
final_df[MZ_COLUMN] = final_df[MZ_COLUMN].round(5)

final_df.to_csv(output_csv, index=False)

print(f"Done. Wrote {len(final_df)} transitions to {output_csv}")