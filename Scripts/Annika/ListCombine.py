import pandas as pd

# Path to your file
file_path = r"Z:\Group Share\Annika\Stellar\EC\Replicates\Listcomparison.xlsx"

# Read Excel with first row as header
df = pd.read_excel(file_path, header=0)

# Clean headers (remove extra spaces)
df.columns = df.columns.str.strip()

# Now your headers are exactly "List 1" and "List 2"
columns_to_combine = ["List 1", "List 2"]

all_values = []
for col in columns_to_combine:
    if col in df.columns:
        col_values = df[col].dropna().astype(str).tolist()
        all_values.extend(col_values)
    else:
        print(f"Warning: column '{col}' not found in Excel sheet.")

# Deduplicate and sort alphabetically
combined = sorted(set(all_values))

# Save combined list
save_path = r"Z:\Group Share\Annika\Stellar\EC\Replicates\combined.xlsx"
out_df = pd.DataFrame(combined, columns=["Combined_List"])
out_df.to_excel(save_path, index=False)

print(f"Combined list saved to {save_path}")
