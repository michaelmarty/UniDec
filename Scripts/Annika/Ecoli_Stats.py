import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
import os

# ----------------------------
# USER INPUTS
# ----------------------------
input_xlsx = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\Normalized_NewColumn_Areas_1,2,3,4,5,6.xlsx"
output_dir = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data"
os.makedirs(output_dir, exist_ok=True)

# Only sheets 1–6
sheets_to_use = ['1', '2', '3', '4', '5', '6']

# Classes to analyze
CLASSES = ['CL', 'PE', 'PG']

# Mapping: controls vs samples (paired)
CONTROL_SHEETS = ['1', '3', '5']
SAMPLE_SHEETS  = ['2', '4', '6']

# ----------------------------
# LOAD AND CLEAN DATA
# ----------------------------
data = {}

for sheet in sheets_to_use:
    df = pd.read_excel(input_xlsx, sheet_name=sheet)

    # Clean column names
    df.columns = df.columns.str.strip().str.replace('\xa0', '', regex=False)

    # Automatically pick first column as lipid name
    lipid_col = df.columns[0]
    df = df.rename(columns={lipid_col: "Lipid"})

    # Extract class
    df["Class"] = df["Lipid"].astype(str).str.split().str[0]

    # Keep only IS-normalized column
    if "IS_normalized_class" not in df.columns:
        raise ValueError(f"IS_normalized_class not found in sheet {sheet}")

    data[sheet] = df[['Lipid', 'Class', 'IS_normalized_class']]

# ----------------------------
# COMBINE DATA LONG FORMAT
# ----------------------------
records = []

for sheet, df in data.items():
    for _, row in df.iterrows():
        records.append({
            "Sheet": sheet,
            "Lipid": row["Lipid"],
            "Class": row["Class"],
            "IS_norm": row["IS_normalized_class"]
        })

long_df = pd.DataFrame(records)

# ----------------------------
# LIPID-LEVEL STATS
# ----------------------------
lipid_stats = []

for lipid, sub in long_df.groupby("Lipid"):
    ctrl_vals = []
    samp_vals = []

    for c, s in zip(CONTROL_SHEETS, SAMPLE_SHEETS):
        ctrl_val = sub.loc[sub["Sheet"] == c, "IS_norm"].values
        samp_val = sub.loc[sub["Sheet"] == s, "IS_norm"].values

        if len(ctrl_val) == 1 and len(samp_val) == 1:
            ctrl_vals.append(ctrl_val[0])
            samp_vals.append(samp_val[0])

    if len(ctrl_vals) < 2:
        continue  # skip if not enough pairs

    ctrl_vals = np.array(ctrl_vals)
    samp_vals = np.array(samp_vals)

    fc = np.mean(samp_vals) / np.mean(ctrl_vals)
    log2fc = np.log2(fc)
    tstat, pval = ttest_rel(samp_vals, ctrl_vals)

    lipid_stats.append({
        "Lipid": lipid,
        "Class": sub["Class"].iloc[0],
        "Ctrl_mean": np.mean(ctrl_vals),
        "Samp_mean": np.mean(samp_vals),
        "Ctrl_std": np.std(ctrl_vals, ddof=1),
        "Samp_std": np.std(samp_vals, ddof=1),
        "FoldChange": fc,
        "Log2FC": log2fc,
        "PValue": pval,
        "NegLog10P": -np.log10(pval)
    })

lipid_stats_df = pd.DataFrame(lipid_stats)
lipid_stats_df.to_csv(os.path.join(output_dir, "lipid_level_stats.csv"), index=False)

# ----------------------------
# LIPID-LEVEL VOLCANO PLOTS
# ----------------------------
# Define colors for each class
class_colors = {
    "CL": "dodgerblue",
    "PE": "green",
    "PG": "orange"
}

plt.figure(figsize=(7,5))

for cls, color in class_colors.items():
    sub = lipid_stats_df[lipid_stats_df["Class"] == cls]
    plt.scatter(sub["Log2FC"], sub["NegLog10P"], label=cls, color=color)

plt.axvline(0, linestyle="--", color='gray')
plt.axhline(-np.log10(0.05), linestyle="--", color='red')
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.title("Volcano Plot - All Classes")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "volcano_all_classes.png"), dpi=300)
plt.close()

# ----------------------------
# CLASS-LEVEL STATS
# ----------------------------
class_stats = []

for cls in CLASSES:
    ctrl_sums = []
    samp_sums = []

    for c, s in zip(CONTROL_SHEETS, SAMPLE_SHEETS):
        ctrl_sum = long_df[(long_df["Sheet"] == c) & (long_df["Class"] == cls)]["IS_norm"].sum()
        samp_sum = long_df[(long_df["Sheet"] == s) & (long_df["Class"] == cls)]["IS_norm"].sum()

        ctrl_sums.append(ctrl_sum)
        samp_sums.append(samp_sum)

    fc = np.mean(samp_sums) / np.mean(ctrl_sums)

    class_stats.append({
        "Class": cls,
        "Ctrl_mean": np.mean(ctrl_sums),
        "Samp_mean": np.mean(samp_sums),
        "Ctrl_std": np.std(ctrl_sums, ddof=1),
        "Samp_std": np.std(samp_sums, ddof=1),
        "FoldChange": fc,
        "Log2FC": np.log2(fc)
    })

class_stats_df = pd.DataFrame(class_stats)
class_stats_df.to_csv(os.path.join(output_dir, "class_level_stats.csv"), index=False)

# ----------------------------
# CLASS-LEVEL FOLD-CHANGE BAR PLOT
# ----------------------------
plt.figure(figsize=(5,4))
plt.bar(class_stats_df["Class"], class_stats_df["Log2FC"], color='salmon')
plt.axhline(0, linestyle="--", color='gray')
plt.ylabel("Log2 Fold Change")
plt.title("Class-Level Fold Changes")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "class_fold_change.png"), dpi=300)
plt.close()

print("✅ Lipid and class stats + volcano/fold-change plots complete!")