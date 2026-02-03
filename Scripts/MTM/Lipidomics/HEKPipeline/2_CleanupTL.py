import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp
from PublicScripts.Lipids.LipidFunctions import rename_overlapping_species

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Final_results"
os.chdir(topdir)

filepath = r"Combined_TransitionList.csv"
outfile = filepath.replace(".csv", "_Clean.csv")

df = pd.read_csv(filepath)

# Remove any where "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

# Drop classes with these names
drop_classes = ["ASM", "EtherSMGDG", "EtherTG", "MGDG", "EtherMGDG", "LNAPS", "LNAPE", "DLCL", "EtherLPG", "MLCL"]
df_filtered = df_filtered[~df_filtered["Molecule List Name"].isin(drop_classes)]

# Only take PC class for testing
# df_filtered = df_filtered[df_filtered["Molecule List Name"] == "PC"]
# # only take [M+CH3COO] adduct
# df_filtered = df_filtered[df_filtered["Precursor Adduct"] == "[M+H]+"]
# # only take name contains PC 16:1_18:1
# df_filtered = df_filtered[df_filtered["Molecule Name"].str.contains("PC 18:1_22:1")]

# Remove any where Precursor Adduct is [M-2H]
df_filtered = df_filtered[df_filtered["Precursor Adduct"] != "[M-2H]"]

# Remove any where Class is PG and Product Mz is close to 227.0326 (common precursor fragment)
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(227.0326 - 0.1, 227.0326 + 0.1)))]
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(171.006 - 0.1, 171.006 + 0.1)))]

# Remove PE fragment near 152.9958
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PE") & (df_filtered["Product Mz"].between(152.9958 - 0.1, 152.9958 + 0.1)))]

df_filtered = lfp.skyline_tail_namer(df_filtered)

# if "Molecule Formula" in df_filtered.columns:
#     df_filtered = df_filtered.drop(columns=["Molecule Formula"])

if "Molecule Formula" in df_filtered.columns:
    # Drop any rows with "[2H]" in the formula
    df_filtered = df_filtered[~df_filtered["Molecule Formula"].str.contains(r"\[2H\]")]

# If charge is 0 set to 1
# df_filtered.loc[df_filtered["Precursor Charge"] == 0, "Precursor Charge"] = 1
# If charge is 0 drop those rows
df_filtered = df_filtered[df_filtered["Precursor Charge"] != 0]

# Rename species with the same name but different RTs
df_filtered = rename_overlapping_species(df_filtered)

# Drop any fragments where Product Name contains +HG
# df_filtered = df_filtered[~df_filtered["Product Name"].str.contains(r"\+HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("NOT")]

df_filtered.to_csv(outfile, index=False)