import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp
from PublicScripts.Lipids.LipidFunctions import rename_overlapping_species

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2"
os.chdir(topdir)

filepath = r"Combined_TransitionList_Cleanest_3.csv"
outfile = filepath.replace("_3.csv", "_4.csv")

df = pd.read_csv(filepath)

# Remove any where "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

# Drop classes with these names
drop_classes = ["ASM", "EtherSMGDG", "EtherTG", "MGDG", "EtherMGDG", "LNAPS", "LNAPE", "DLCL", "EtherLPG", "MLCL",
                "AHexCer", "HexCer_NDS", "EtherPS", "OxPE", "LPA", "LPI", "LPS", "EtherDG", "CoQ", "DGDG", "NATau", "PE_Cer", "HBMP", "EtherLPC", "EtherPC", "PEtOH"]
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
              "TG 15:0_16:0_18:1", "TG 16:0_18:0_32:0", "Cer 30:1;2|Cer 15:1;2/15:0", "Cer 32:1;2|Cer 18:1;2/14:0", "Cer 34:2;2|Cer 18:2;2/16:0", "Cer 35:1;2|Cer 18:1;2/17:0",
              "Cer 36:2;2|Cer 18:1;2/18:1", "Cer 36:2;2|Cer 18:2;2/18:0", "Cer 38:2;2|Cer 18:2;2/20:0", "Cer 39:3;2|Cer 18:1;2/21:2","Cer 40:2;2|Cer 18:2;2/22:0",
              "Cer 43:1;2|Cer 18:1;2/25:0","Cer 44:0;2|Cer 18:0;2/26:0" "Hex2Cer 34:1;2", "Hex2Cer 38:1;2|Hex2Cer 17:0;2/21:1", "Hex2Cer 40:2;2|Hex2Cer 18:1;2/22:1", "Hex2Cer 42:2;2",
              "HexCer 32:1;2|HexCer 18:1;2/14:0_A", "HexCer 32:1;2|HexCer 18:1;2/14:0_B", "HexCer 32:2;2|HexCer 18:1;2/14:1", "HexCer 40:2;2|HexCer 18:2;2/22:0", "HexCer 42:2;2|HexCer 18:2;2/24:0_A", "HexCer 42:2;2|HexCer 18:2;2/24:0_B",
              "LPC 15:0", "LPC 16:2", "LPC 17:2", "LPC 18:2", "LPC 19:1", "LPC 19:2", "LPC 20:2", "LPC 20:3", "LPC 22:2", "LPC 22:3", "LPC 24:2", "LPC 28:2", "LPE 14:0", "LPE 15:0", "LPE 16:2", "LPE 18:2", "LPE 18:3", "LPE 19:0", "LPE 19:1", "LPE 19:2", "LPE 20:2", "LPE 20:3",
              "LPE 20:4", "LPE 21:1", "LPE 22:1", "LPE 22:2", "LPE 22:3", "LPE 24:0", "LPE 24:1", "LPE 24:2", "LPE 24:3", "LPE 28:1", "LPE O_20:5", "LPE O_22:4", "LPE O_22:6", "LPG 15:0", "LPG 16:1", "LPG 20:1", "LPG  20:3_A", "LPG 20:3_B", "PA 38:2", "PC 20:0", "PC 29:2_A", "PC 29:2_B",
              "PC 30:2", "PC 31:1", "PC 33:0", "PC 33:1", "PC 33:2", "PC 33:3", "PC 33:6", "PC 34:4", "PC 35:0", "PC 35:1", "PC 35:2", "PC 35:3", "PC 35:6", "PC 36:5", "PC 36:7", "PC 37:2", "PC 37:5", "PC 37:9", "PC 38:1", "PC 39:9", "PC 40:0", "PC 40:6", "PC 41:0", "PC 41:7", "PC 41:9",
              "PC 42:1", "PC 43:2", "PC 44:1", "PC 47:2", "PC 49:1", "PC 50:1", "PC 50:3", "PC 51:7", "PC 52:4", "PC 54:1", "PC 54:4", "PC 55:8", "PC 56:2", "PC O_17:1", "PC O_22:0","PC O_29:1", "PC O_30:0", "PC O_31:0", "PC O_31:1", "PC O_31:2",
              "PC O_32:0", "PC O_32:1", "PC O_33:0_A", "PC O_33:0_B", "PC O_33:1", "PC O_33:2_A", "PC O_33:2_B", "PC O_33:5", "PC O_34:0", "PC O_34:5", "PC O_34:6", "PC O_34:7", "PC O_34:9", "PC O_35:0", "PC O_36:10",  "PC O_35:7", "PC O_36:5", "PC O_36:9", "PC O_37:1", "PC O_38:0", "PC O_38:1","PC O_38;2", "PC O_39:1","PC O_39:1",
              "PC O_39:2", "PC O_39:3", "PC O_40:0", "PC O_41:3", "PC O_42:7", "PC O_40:1", "PC O_40:2", "PC O_41:1", "PC O_42:1", "PC O_42:2", "PC O_44:2", "PC O_45:1", "PC O_47:1", "PC O_49:1", "PC O_49:2", "PC O_51:1", "PC O_51:2", "PC O_51:3", "PC O_53:1", "PC O_53:2", "PC O_53:3",
              "PC O_54:1","PC O_54:3", "PC O_55:2", "PC O_55:3", "Cer 44:3;4|Cer 28:2;3(FA 16:0)", "Cer 44:4;4|Cer 26:2;3(FA 18:1)", "Cer 48:4;4|Cer 28:2;3(FA 20:1)", "CL 64:1|CL 14:0_16:0_16:0_18:1", "LPC 24:1", "PC O_16:0", "PC O_19:3", "PC O_24:0",
              "PC O_26:1", "PC O_32:5", "PC O_32:6", "PC O_33:5", "PE O_18:0", "SM 18:0;O2", "SM 40:1;O2", "SM 17:1;2", "SM 19:0;2", "SM 21:1;2", "SM 23:0;2", "SM 34:0;2", "SM 36:1;2|SM 18:1;2/18:0", "PE 36:1", "PE 40:2", "PE 42:2","PE 48:1",
              "PE 48:2", "PE 50:1", "PE 54:3", "PI 38:2", "PI 38:4", "PI 40:2", "PI 40:3", "PS 15:0", "PS 17:1", "PS 32:1", "PS 40:2", "TG 42:1|TG 12:0_12:0_18:1", "TG 44:0|TG 14:0_14:0_16:0", "PC 26:0", "PC 28:1", "PC 29:0", "PC 30:1", "PC 31:0",
              "SM 42:1;2", "SM 42:2;2", "SM 42:3;2", "SM 43:7;2", "SM 44:2;2", "SM 47:0;2", "SM 49:1;2", "SM 49:2;2", "SM 51:1;2", "SM 53:1;2", "SM 53:2;2", "SM 53:3;2", "SM 54:2;2", "SM 55:0;2", "SM 55:1;2", "SM 55:2;2", "SM 55:3;2", "SM 56:0;2",
              "SM 57:3;2", "SM 60:1;2", "SM 31:0;2", "SM 33:0;2", "SM 33:1;2", "SM 33:3;2", "SM 35:0;2", "SM 40:5;2|SM 18:2;2/22:3", "SM 43:2;2|SM 17:0;2/26:2", "SM 52:2;2|SM 18:1;2/34:1", "SM 57:1;2|SM 18:1;2/39:0", "SM 58:2;2|SM 18:0;2/40:2", "SM 59:1;2|SM 18:1;2/41:0",
              "SM 61:1;2|SM 17:0;2/44:1", "SM 61:2;2|SM 17:0;2/44:2", "SM 61:3;2|SM 29:3;2/32:0", "SM 37:5;2", "SM 37:6;2", "SM 39:3;2", "SM 39:6;2", "SM 35:2;2", "SM 37:1;2", "SM 39:1;2",
]
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

# if "Molecule Formula" in df_filtered.columns:
#     # Drop any rows with "[2H]" in the formula
#     df_filtered = df_filtered[~df_filtered["Molecule Formula"].str.contains(r"\[2H\]")]

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