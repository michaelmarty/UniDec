import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp
from PublicScripts.Lipids.LipidFunctions import rename_overlapping_species


filepath = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Full Library Master Lists\Positive_updated_PeakMaster_7_2026_01_23_14_30_21_TransitionList.csv"
outfile = filepath.replace(".csv", "_NoPrecursorAssigned.csv")

df = pd.read_csv(filepath)

# Remove anywhere "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

# Only take PC class for now
df_filtered = df_filtered[df_filtered["Molecule List Name"] == "PC"]

df_filtered = rename_overlapping_species(df_filtered)

# Remove anywhere Precursor Adduct is [M-2H]
# df_filtered = df_filtered[df_filtered["Precursor Adduct"] != "[M-2H]"]

# Remove anywhere Class is PG and Product Mz is close to 227.0326 (common precursor fragment)
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(227.0326 - 0.1, 227.0326 + 0.1)))]
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(171.006 - 0.1, 171.006 + 0.1)))]

# Remove PE fragment near 152.9958
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PE") & (df_filtered["Product Mz"].between(152.9958 - 0.1, 152.9958 + 0.1)))]

df_filtered = lfp.skyline_tail_namer(df_filtered)

# Drop any fragments where Product Name contains +HG
# df_filtered = df_filtered[~df_filtered["Product Name"].str.contains(r"\+HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("NOT")]

# Drop molecular formula column if it exists
if "Molecule Formula" in df_filtered.columns:
    df_filtered = df_filtered.drop(columns=["Molecule Formula"])

# If charge is 0 set to 1
# df_filtered.loc[df_filtered["Precursor Charge"] == 0, "Precursor Charge"] = 1
# If charge is 0 drop those rows
df_filtered = df_filtered[df_filtered["Precursor Charge"] != 0]



df_filtered.to_csv(outfile, index=False)