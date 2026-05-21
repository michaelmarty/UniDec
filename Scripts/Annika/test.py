import pandas as pd

# === Update these paths ===
iso_file = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData\HEKExtr_Fixed_Clean16_isolist.csv"  # has m/z, z, Polarity, RT Window
skyline_rt_file = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData\Fixed_Windows_IsoList.csv"      # Skyline export: Compound + RT
output_file = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData\Final_Isolation_List.csv"

compound_col = "Compound"
rt_col = "Retention Time (min)"  # change if Skyline uses a different header

# === Read files ===
iso_df = pd.read_csv(iso_file)
rt_df = pd.read_csv(skyline_rt_file)

# Clean column names
iso_df.columns = iso_df.columns.str.strip()
rt_df.columns = rt_df.columns.str.strip()

# Clean compound names
iso_df[compound_col] = iso_df[compound_col].astype(str).str.strip()
rt_df[compound_col] = rt_df[compound_col].astype(str).str.strip()

# If Skyline has duplicates per compound, keep first
rt_df_unique = rt_df.drop_duplicates(subset=[compound_col], keep="first")

# === Merge: keep all iso list rows, add/replace RT from Skyline ===
# First, rename Skyline RT column to avoid confusion
rt_df_unique = rt_df_unique[[compound_col, rt_col]].rename(columns={rt_col: "RT_from_Skyline"})

merged = iso_df.merge(rt_df_unique, on=compound_col, how="left")

# If you want to REPLACE the iso RT with Skyline RT where available:
if rt_col in merged.columns:
    merged[rt_col] = merged["RT_from_Skyline"].combine_first(merged[rt_col])
else:
    # If iso file didn't have RT column, just create it
    merged[rt_col] = merged["RT_from_Skyline"]

# Drop helper column
merged = merged.drop(columns=["RT_from_Skyline"])

# === Save final isolation list ===
merged.to_csv(output_file, index=False)
print("Saved final isolation list to:", output_file)