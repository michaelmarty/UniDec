import pandas as pd

# Input / output file paths
input_xls = r"Z:\Group Share\Annika\UltimateSplash_TransList.xls"
output_csv = r"Z:\Group Share\Annika\Fixed_UltimateSplash_TransList.csv"
output_xls = r"Z:\Group Share\Annika\Fixed_UltimateSplash_TransList.xlsx"  # Excel output

# Read Excel as text to prevent formulas being evaluated
df = pd.read_excel(input_xls, dtype=str)  # ensures all values are read as literal text

# Add "IS " to start of each molecule name
df["PrecursorName"] = "IS " + df["PrecursorName"].astype(str)

# Save as CSV (Excel will treat all values as text if opened via Import)
df.to_csv(output_csv, index=False)

# Save as XLSX (formulas replaced with literal text)
df.to_excel(output_xls, index=False, engine='openpyxl')

print("Done! Saved CSV to:", output_csv)
print("Done! Saved XLSX to:", output_xls)
