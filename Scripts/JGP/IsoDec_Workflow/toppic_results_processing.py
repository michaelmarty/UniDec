import pandas as pd
import numpy as np

df = pd.read_csv("D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\20170105_L_MaD_ColC4_Ecoli20161116-12_01_ms2_toppic_prsm.tsv", sep="\t", skiprows=49)

total_prsms = 0
filtered_prsms = 0

raw_matched_frags = []
filtered_matched_frags = []

raw_evals = []
filtered_evals = []

for i, row in df.iterrows():
    total_prsms += 1
    if row["#unexpected modifications"] == 0:
        filtered_prsms += 1
        if not np.isnan(row["#matched fragment ions"]):
            filtered_matched_frags.append(row["#matched fragment ions"])
        if not np.isnan(row["E-value"]):
            filtered_evals.append(row["E-value"])
    if not np.isnan(row["#matched fragment ions"]):
        raw_matched_frags.append(row["#matched fragment ions"])
    if not np.isnan(row["E-value"]):
        raw_evals.append(row["E-value"])

print("Total PRSMs: ", total_prsms)
print("Filtered PRSMs: ", filtered_prsms)
print("Average matched fragments (raw): ", np.mean(raw_matched_frags))
print("Average matched fragments (filtered): ", np.mean(filtered_matched_frags))
print("Average E-value (raw): ", np.mean(raw_evals))
print("Average E-value (filtered): ", np.mean(filtered_evals))
