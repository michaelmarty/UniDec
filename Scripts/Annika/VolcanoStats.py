import pandas as pd
import numpy as np
import os

filename = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\20251216_MoreLipidDiscs_FTs\DDM Comparisons\1%_DDM_Comparison_Stats.xlsx"

df = pd.read_excel(filename)
print(df)

c1 = ["1% 1", "1% 2", "1% 3"]
c2 = ["D1", "D2", "D3"]

# Check normalization
for c in c1:
    print(np.sum(df[c]))
for c in c2:
    print(np.sum(df[c]))


# Run a paired t-test between the two conditions for each metric
from scipy.stats import ttest_rel
results = []
for index, row in df.iterrows():
    vals1 = row[c1].to_numpy(dtype=float)
    vals2 = row[c2].to_numpy(dtype=float)
    tstat, pval = ttest_rel(vals1, vals2)
    foldchanges = vals2 / vals1
    meanfoldchange = np.mean(foldchanges)
    results.append({
        "Fold Change (DDM/1% DDM)": meanfoldchange,
        "T-statistic": tstat,
        "P-value": pval
    })

# Perform B-H correction on p-values
pvals = [res["P-value"] for res in results]
from statsmodels.stats.multitest import multipletests
reject, pvals_corrected, _, _ = multipletests(pvals, method="fdr_bh")
for i, res in enumerate(results):
    res["Adjusted P-value"] = pvals_corrected[i]
    res["Significant"] = reject[i]

# Add results to df
results_df = pd.DataFrame(results)
final_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)
# Save to new Excel file
output_filename = filename.replace(".xlsx", "_Stats.xlsx")
final_df.to_excel(output_filename, index=False)
