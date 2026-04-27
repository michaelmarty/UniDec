import pandas as pd
import csv
# === SETTINGS ===
input_xls = r"Z:\Group Share\Annika\UltimateSplash_TransList.xls"
output_csv = r"Z:\Group Share\Annika\New_UltimateSplash_TransList_Clean.csv"
column_name = "ProductName"  # change this


# Read Excel as text (prevents numeric evaluation)
df = pd.read_excel(input_xls, sheet_name=0, dtype=str)

# Replace any cell starting with '=' with its literal text (remove the '=')
df = df.map(lambda x: x[1:] if isinstance(x, str) and x.startswith('=') else x)

# Save to CSV
df.to_csv(output_csv, index=False)

print(f"Saved CSV to {output_csv} without formulas")