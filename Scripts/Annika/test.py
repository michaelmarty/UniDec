import pandas as pd

# ----------------------------
# User settings
# ----------------------------
input_excel = r"Z:\Group Share\Annika\Stellar\HEK\TargetedHEK_StellarList.xlsx"  # your Excel file
sheet_all = "TargetedHEK_StellarList"      # sheet with all lipids
sheet_pos = "Pos"       # sheet with positive lipids
output_excel = r"Z:\Group Share\Annika\Stellar\HEK\Lipids_with_Overlap.xlsx"
key_column = "Compound"  # column to match lipids
# ----------------------------

# Load the sheets
all_df = pd.read_excel(input_excel, sheet_name=sheet_all)
pos_df = pd.read_excel(input_excel, sheet_name=sheet_pos)

# Ensure key column exists
for df, name in [(all_df, sheet_all), (pos_df, sheet_pos)]:
    if key_column not in df.columns:
        raise ValueError(f"Column '{key_column}' not found in sheet '{name}'")

# Find overlapping lipids (keep all columns from AllLipids)
overlap_df = all_df[all_df[key_column].isin(pos_df[key_column])]

print(f"Found {len(overlap_df)} overlapping lipids.")

# Write back to Excel with a new sheet (keep original sheets)
with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    all_df.to_excel(writer, sheet_name=sheet_all, index=False)
    pos_df.to_excel(writer, sheet_name=sheet_pos, index=False)
    overlap_df.to_excel(writer, sheet_name="Overlap", index=False)

print(f"Overlap sheet added to {output_excel}")
