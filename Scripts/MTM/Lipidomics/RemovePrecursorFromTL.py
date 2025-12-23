import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp

filepath = r"Z:\Group Share\Dash\Ecoli_TransList_FINAL.csv"
outfile = filepath.replace(".csv", "_NoPrecursorAssigned.csv")

df = pd.read_csv(filepath)

# Remove any where "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

df_filtered = lfp.skyline_tail_namer(df_filtered)

df_filtered.to_csv(outfile, index=False)