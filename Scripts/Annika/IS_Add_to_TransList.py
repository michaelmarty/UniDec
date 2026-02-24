import pandas as pd

# Input / output file paths
input_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\UltimateSPLASH_TGs_TransList_Clean.csv"
output_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Fixed_UltimateSPLASH_TGs_TransList_Clean.csv"


# Read Excel as text to prevent formulas being evaluated
df = pd.read_csv(input_csv, dtype=str)  # ensures all values are read as literal text

# Add "IS " to start of each molecule name
df["Molecule Name"] = "IS " + df["Molecule Name"].astype(str)

# Save as CSV (Excel will treat all values as text if opened via Import)
df.to_csv(output_csv, index=False)

print("Done! Saved CSV to:", output_csv)
