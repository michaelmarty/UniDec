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
drop_classes = ["ASM", "EtherSMGDG", "EtherTG", "MGDG", "EtherMGDG", "LNAPS", "LNAPE", "DLCL", "EtherLPG", "MLCL",
                "AHexCer", "HexCer_NDS", "EtherPS", "OxPE", "LPA", "LPI"]
df_filtered = df_filtered[~df_filtered["Molecule List Name"].isin(drop_classes)]

# Drop TG, LPC, PC, SM with M+Na adducts
df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "TG") & (df_filtered["Precursor Adduct"] == "[M+Na]+"))]
df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "SM") & (df_filtered["Precursor Adduct"] == "[M+Na]+"))]
df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PC") & (df_filtered["Precursor Adduct"] == "[M+Na]+"))]
df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "LPC") & (df_filtered["Precursor Adduct"] == "[M+Na]+"))]

# On LPC, strip /0:0 from names
df_filtered.loc[df_filtered["Molecule List Name"] == "LPC", "Molecule Name"] = df_filtered.loc[df_filtered["Molecule List Name"] == "LPC", "Molecule Name"].str.replace(r"/0:0", "", regex=True)


# Group Cer classes together for Cer_NDS, Cer_EOS, Cer_NS, Cer_EODS
cer_classes = ["Cer_NDS", "Cer_EOS", "Cer_NS", "Cer_EODS"]
df_filtered.loc[df_filtered["Molecule List Name"].isin(cer_classes), "Molecule List Name"] = "Cer"

# Only take PC class for testing
# df_filtered = df_filtered[df_filtered["Molecule List Name"] == "Cer_NS"]

# Drop specific problems
drop_names = ["CL 14:0_18:1_16:1_16:1", r"Cer 19:0;2O/43:0", r"Cer 18:1;O2/40:2", r"Cer 18:1;O2/40:3", r"Cer 18:1;O2/42:3",
              "Cer 56", "Cer 58", "Cer 60", "Cer 62", "Cer 52", "Cer 54", "Cer 51", "PC 85:4", "PC 57:5",
              "TG 15:0_16:0_18:1", "TG 16:0_18:0_32:0"]
# Drop if Molecule Name contains any of the drop names
for name in drop_names:
    df_filtered = df_filtered[~df_filtered["Molecule Name"].str.contains(name)]
# print(df_filtered.to_string())
# exit()

# Drop specific transitions in specific lipids
drop_frags_dict = {"PS 16:0_16:1":[253.2166], "PS 16:0_18:1" : [281.2479]}
for lipid, frags in drop_frags_dict.items():
    for frag in frags:
        df_filtered = df_filtered[~((df_filtered["Molecule Name"].str.contains(lipid)) & (df_filtered["Product Mz"].between(frag - 0.1, frag + 0.1)))]

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
df_filtered = rename_overlapping_species(df_filtered, rttol=1)

df_filtered = lfp.skyline_tail_namer(df_filtered)

# Drop any fragments where Product Name contains +HG
# df_filtered = df_filtered[~df_filtered["Product Name"].str.contains(r"\+HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("NOT")]

df_filtered.to_csv(outfile, index=False)

# PROBLEM LIST
#

"""
Notes:
LPC isomers

"""