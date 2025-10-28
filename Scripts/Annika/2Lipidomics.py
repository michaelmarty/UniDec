import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load your data
# -----------------------------
folder = r"Z:\Group Share\Annika\Stellar\EC\Replicates"
file_path = f"{folder}\\AqpZ_Memb_Ext_Disc_Merged.csv"

# Read CSV, header is row 4 (zero-indexed 3)
df = pd.read_csv(file_path, header=None)


# Clean column names: strip spaces and quotes
df.columns = ["Lipid Name", "ME1", "ME2", "ME3", "D1", "D2", "D3"]
for col in ["ME1","ME2","ME3","D1","D2","D3"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# -----------------------------
# 2. Identify columns
# -----------------------------
feature_col = "Lipid Name"


# Replace these with your actual replicate column names from the CSV
ecoli_cols = ["ME1", "ME2", "ME3"]      # membrane extract replicates
membrane_cols = ["D1", "D2", "D3"]      # disc replicates

# Ensure replicate columns exist
for col in ecoli_cols + membrane_cols:
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in CSV. Available columns: {df.columns.tolist()}")

# Convert replicate columns to numeric, remove quotes if present
for col in ecoli_cols + membrane_cols:
    df[col] = pd.to_numeric(df[col].astype(str).str.replace(r"[\"']", "", regex=True), errors="coerce")

# -----------------------------
# 3. Normalize (median scaling per column)
# -----------------------------
for col in ecoli_cols + membrane_cols:
    df[col] = df[col] / df[col].median()

# -----------------------------
# 4. Fold change + t-test
# -----------------------------
results = []

for _, row in df.iterrows():
    ecoli_vals = row[ecoli_cols].values
    mem_vals = row[membrane_cols].values

    # Safety check for empty or zero slices
    if len(ecoli_vals) == 0 or len(mem_vals) == 0 or np.all(np.isnan(ecoli_vals)) or np.all(np.isnan(mem_vals)):
        log2fc = np.nan
        pval = np.nan
    else:
        mean_ecoli = np.mean(ecoli_vals)
        mean_mem = np.mean(mem_vals)
        log2fc = np.log2(mean_mem / mean_ecoli) if mean_ecoli > 0 else np.nan
        tstat, pval = stats.ttest_ind(mem_vals, ecoli_vals, equal_var=False)

    results.append({
        "Feature": row[feature_col],
        "Class": row[class_col] if class_col in df.columns else "",
        "log2FC": log2fc,
        "pval": pval
    })

res_df = pd.DataFrame(results)
res_df["-log10p"] = -np.log10(res_df["pval"])

# -----------------------------
# 5. Volcano plot
# -----------------------------
plt.figure(figsize=(8, 6))
plt.scatter(res_df["log2FC"], res_df["-log10p"], alpha=0.7)

# Threshold lines
plt.axvline(x=1, color="red", linestyle="--")
plt.axvline(x=-1, color="red", linestyle="--")
plt.axhline(y=-np.log10(0.05), color="blue", linestyle="--")

plt.xlabel("log2 Fold Change (Disc / Membrane Extract)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot: Disc vs Membrane Extract")
plt.show()

# -----------------------------
# 6. Class-level fold change
# -----------------------------
if class_col in df.columns:
    class_fc = (
        df.groupby(class_col)[ecoli_cols + membrane_cols]
        .mean()
        .apply(lambda x: np.log2(np.mean(x[membrane_cols]) / np.mean(x[ecoli_cols])), axis=1)
    )

    plt.figure(figsize=(8, 5))
    class_fc.sort_values().plot(kind="bar", color="gray")
    plt.axhline(y=0, color="black", linestyle="--")
    plt.ylabel("log2 Fold Change (Disc / Membrane Extract)")
    plt.title("Class-level Fold Change")
    plt.tight_layout()
    plt.show()
