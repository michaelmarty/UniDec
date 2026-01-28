import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from pathlib import Path

# -----------------------------
# Load CSV
# -----------------------------
file_path = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\Analysis\NEWSamples1,2,3,4,5,6_ANALYSIS.csv"
df = pd.read_csv(file_path)

# -----------------------------
# Define the groups (replace with your column names)
# -----------------------------
groupA = ["1", "3", "5"]  # columns for condition A
groupB = ["2", "4", "6"]  # columns for condition B

# -----------------------------
# Clean & convert to numeric
# -----------------------------
for col in groupA + groupB:
    df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '').str.strip(), errors='coerce')

# -----------------------------
# Compute t-test p-values row by row
# -----------------------------
pvals = []
for _, row in df.iterrows():
    a = row[groupA].dropna().astype(float).values
    b = row[groupB].dropna().astype(float).values

    if len(a) > 0 and len(b) > 0:
        p = ttest_ind(a, b, equal_var=False).pvalue
    else:
        p = np.nan

    pvals.append(p)

df["PValue"] = pvals

# -----------------------------
# Save results to a CSV you control
# -----------------------------
# Change this path to wherever you want the CSV to appear
save_folder = Path(r"Downloads")  # for example
save_folder.mkdir(parents=True, exist_ok=True)  # ensure folder exists
output_file = save_folder / "PValues_Lipid_Test.csv"

df.to_csv(output_file, index=False)
print(f"Done! Saved p-values to: {output_file}")