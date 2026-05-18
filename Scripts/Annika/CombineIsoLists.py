import pandas as pd

# Input CSV files
file1 = r"Z:\Group Share\Annika\Stellar\Brain Library\Test_BrainAXLib_BrainIsoList.csv"
file2 = r"Z:\Group Share\Annika\Stellar\Brain Library\Test_HEKLib_BrainIsoList.csv"

# Output file
output_file = r"Z:\Group Share\Annika\Stellar\Brain Library\CombinedIsoList_Brain.csv"

# Read CSVs
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Combine the two dataframes
combined = pd.concat([df1, df2], ignore_index=True)

# Remove duplicate compounds
# Keeps the first occurrence of each compound
deduplicated = combined.drop_duplicates(subset=["Compound"])

# Optional: sort by compound name
deduplicated = deduplicated.sort_values(by="Compound")

# Save to new CSV
deduplicated.to_csv(output_file, index=False)

print(f"Deduplicated CSV saved as: {output_file}")