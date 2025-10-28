import pandas as pd

# ----------------------------
# User settings
# ----------------------------
file_path = r"Z:\Group Share\Annika\Stellar\HEK\LipidCreatorTransitionList.csv"
output_file = r"Z:\Group Share\Annika\Stellar\HEK\precursor_only.csv"
# ----------------------------

# Load the transition list (CSV)
df = pd.read_csv(file_path)  # if your file is CSV
# df = pd.read_excel(file_path, header=0)  # use this if it's XLSX

# Check for common column names
# We'll look for these typical columns: PrecursorName, PrecursorMz, PrecursorAdduct, PrecursorCharge
required_columns = ['PrecursorName', 'PrecursorMz']
for col in required_columns:
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in input file. Columns found: {df.columns.tolist()}")

# Keep only the unique precursors
precursor_df = df.drop_duplicates(subset=['PrecursorName', 'PrecursorMz'])

# Optional: keep adduct and charge if they exist
for optional_col in ['PrecursorAdduct', 'PrecursorCharge']:
    if optional_col in df.columns:
        precursor_df = precursor_df[[*precursor_df.columns, optional_col]].drop_duplicates()

# Keep only the main columns in a clean order
output_columns = ['PrecursorName', 'PrecursorMz']
for col in ['PrecursorAdduct', 'PrecursorCharge']:
    if col in precursor_df.columns:
        output_columns.append(col)
precursor_df = precursor_df[output_columns]

# Save to CSV
precursor_df.to_csv(output_file, index=False)

print(f"Precursor-only list saved to '{output_file}' with {len(precursor_df)} unique precursors.")
