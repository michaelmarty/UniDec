import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp
from PublicScripts.Lipids.LipidFunctions import rename_overlapping_species

topdir = r"Z:\Group Share\Annika\IQ-X\AcquireX and Neo\AcquireX"
os.chdir(topdir)

filepath = r"Brain_PosNeg_TransList.csv"
outfile = filepath.replace("_TransList.csv", "_TransList_Clean1.csv")

df = pd.read_csv(filepath)

# Remove any where "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

# Drop classes with these names
drop_classes = ["ASM", "EtherSMGDG", "EtherTG", "MGDG", "EtherMGDG", "LNAPS", "LNAPE", "DLCL", "EtherLPG", "MLCL",
                "AHexCer", "HexCer_NDS", "EtherPS", "OxPE", "EtherDG", "DGDG", "NATau", "PE_Cer", "HBMP", "PEtOH", "DGGA", "ADGGA", "Cer_AS", "EtherLPE",
                "HexCer_HDS", "HexCer_AP", "HexCer_EOS", "BRSE", "Sph", "SHex", "EtherPG", "EtherPI", "PI_Cer", "PMeOH", "Cer_ADS"

]
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
              "Cer 43:1;2|Cer 18:1;2/25:0","Cer 44:0;2|Cer 18:0;2/26:0", "HexCer 32:1;2|HexCer 18:1;2/14:0_A", "HexCer 32:1;2|HexCer 18:1;2/14:0_B", "HexCer 32:2;2|HexCer 18:1;2/14:1",
              "LPC 15:0", "LPC 17:2", "LPC 19:1", "LPC 19:2", "LPC 24:2", "LPC 28:2", "LPE 14:0", "LPE 15:0", "LPE 16:2", "LPE 18:2", "LPE 18:3", "LPE 19:0", "LPE 19:1", "LPE 19:2", "LPE 20:2", "LPE 20:3",
              "LPE 20:4", "LPE 21:1", "LPE 22:1", "LPE 22:2", "LPE 22:3", "LPE 24:0", "LPE 24:1", "LPE 24:2", "LPE 24:3", "LPE 28:1", "LPE O_20:5", "LPE O_22:4", "LPE O_22:6", "LPG 15:0", "LPG 16:1", "LPG 20:1", "LPG  20:3_A", "LPG 20:3_B", "PA 38:2", "PC 20:0", "PC 29:2_A", "PC 29:2_B",
                "PC 33:6", "PC 35:6",  "PC 36:7", "PC 37:9", "PC 39:9", "PC 41:7", "PC 41:9", "PC 43:2", "PC 44:1", "PC 47:2", "PC 49:1", "PC 50:1", "PC 50:3", "PC 51:7", "PC 52:4", "PC 54:1", "PC 54:4",
              "PC 55:8", "PC 56:2", "PC O_29:1", "PC O_33:5", "PC O_34:5", "PC O_34:6", "PC O_34:7", "PC O_34:9", "PC O_36:10",  "PC O_35:7", "PC O_36:9", "PC O_37:1", "PC O_38:0", "PC O_38:1","PC O_38;2", "PC O_39:1","PC O_39:1",
              "PC O_39:2", "PC O_39:3", "PC O_40:0", "PC O_41:3", "PC O_42:7", "PC O_40:1", "PC O_40:2", "PC O_41:1", "PC O_42:1", "PC O_42:2", "PC O_44:2", "PC O_45:1", "PC O_47:1", "PC O_49:1", "PC O_49:2", "PC O_51:1", "PC O_51:2", "PC O_51:3", "PC O_53:1", "PC O_53:2", "PC O_53:3",
              "PC O_54:1","PC O_54:3", "PC O_55:2", "PC O_55:3", "Cer 44:3;4|Cer 28:2;3(FA 16:0)", "Cer 44:4;4|Cer 26:2;3(FA 18:1)", "Cer 48:4;4|Cer 28:2;3(FA 20:1)", "CL 64:1|CL 14:0_16:0_16:0_18:1", "PC O_19:3", "PC O_24:0",
              "PC O_26:1", "PC O_32:5", "PC O_32:6", "PC O_33:5", "SM 18:0;O2", "SM 17:1;2", "SM 19:0;2", "SM 21:1;2", "SM 23:0;2", "PE 48:1",
              "PE 48:2", "PE 50:1", "PE 54:3", "PS 15:0", "PS 17:1", "PS 32:1", "PS 40:2", "TG 42:1|TG 12:0_12:0_18:1", "TG 44:0|TG 14:0_14:0_16:0", "PC 26:0", "PC 28:1", "PC 29:0",
               "SM 43:7;2", "SM 47:0;2", "SM 49:1;2", "SM 49:2;2", "SM 51:1;2", "SM 53:1;2", "SM 53:2;2", "SM 53:3;2", "SM 54:2;2", "SM 55:0;2", "SM 55:1;2", "SM 55:2;2", "SM 55:3;2", "SM 56:0;2",
              "SM 57:3;2", "SM 60:1;2", "SM 31:0;2", "SM 33:0;2", "SM 33:1;2", "SM 33:3;2", "SM 35:0;2", "SM 43:2;2|SM 17:0;2/26:2", "SM 52:2;2|SM 18:1;2/34:1", "SM 57:1;2|SM 18:1;2/39:0", "SM 58:2;2|SM 18:0;2/40:2", "SM 59:1;2|SM 18:1;2/41:0",
              "SM 61:1;2|SM 17:0;2/44:1", "SM 61:2;2|SM 17:0;2/44:2", "SM 61:3;2|SM 29:3;2/32:0", "SM 37:5;2", "SM 37:6;2", "SM 39:3;2", "SM 39:6;2", "SM 35:2;2", "SM 37:1;2", "SM 39:1;2", "PC O_38:0", "PC O-40:0|PC O-24:0_16:0",
             "Cer 53:3;2|Cer 38:2;2/15:1", "Cer 55:3;4|Cer 40:1;3(FA 15:1)", "Cer 59:3;4|Cer 42:1;3(FA 17:1)", "Cer 59:4;4|Cer 41:3;3(FA 18:0)", "Hex2Cer 55:12;2", "LPC 16:4", "LPC 18:5", "LPC 22:6", "LPC 28:7", "LPE 18:5", "LPE 20:0", "LPE 22:0",
              "LPS 18:3", "LPS 22:1", "LPS 24:0", "LPS 24:1", "PA 28:0|PA 14:0_14:0", "PC 36:9", "PC 37:7", "PC 38:7", "PC 38:8", "PC 39:7", "PC 40:8", "PC 40:9", "PC 41:8", "PC 42:9",
              "PC 43:7", "PC 43:8", "PC O_30:7", "PC O_32:6", "PC O_32:7", "PC O_32:8", "PC O_33:5", "PC O_33:6", "PC O_33:7", "PC O_34:7", "PC O_34:8", "PC O_34:9", "PC O_35:7", "PC O_36:6|PC O_18:5_18:1", "PC O_36:7",
              "PC O_36:8", "PC O_36:9", "PC O_37:6", "PC O_38:10", "PC O_38:8", "PE O_38:9", "PC O_40:10", "PC O_40:11", "PC O_40:8", "PC O_40:9", "PC O_41:8", "PC O_42:10", "PC O_42:9", "PC O_43:7", "PC O_44:10", "PC O_45:9",
              "PE 29:0|PE 10:0_19:0", "PE 32:6", "PE 35:7", "PE 40:8", "PE 40:9|PE 22:4_18:5", "PE O_32:8", "PE O_36:6|PE O_18:5_18:1",
              "PG 32:5", "PG 34:5", "PG 35:5|PG 17:0_18:5", "PG 38:7", "PG 40:8", "PI 38:8", "PI 40:9", "SM 20:0;2", "SM 21:0;2", "SM 25:5;2", "SM 26:6;2", "SM 27:3;2", "SM 27:6;2", "SM 45:10;2", "SM 45:11;2",
              "SM 46:7;2|SM 18:1;2/28:6", "SM 46:8;2", "SM 46:9;2", "SM 47:2;2", "SM 47:7;2", "SM 47:8;2", "SM 47:9;2", "SM 48:7;2|SM 18:1;2/30:6", "SM 48:8;2", "SM 50:1;2|SM 18:1;2/32:0", "SM 50:2;2", "SM 51:0;2", "SM 51:3;2", "SM 52:1;2",
              "SM 53:0;2", "SM 54:0;2", "SM 54:1;2|SM 18:1;2/36:0", "SM 54:5;2|SM 16:3;2/38:2", "SM 58:7;2|SM 36:2;2/22:5", "SM 60:2;2|SM 18:1;2/42:1", "SM 60:7;2", "SM 60:7;2|SM 34:2;2/26:5", "SM 60:8;2|SM 34:3;2/26:5", "SM 62:8;2|SM 24:1;2/38:7",
              "SM 62:9;2", "SM 35:6;2", "SM 36:5;2", "SM 36:6;2", "SM 36:7;2",  "SM 38:6;2", "SM 38:7;2", "SM 39:7;2", "SM 40:6;2", "SM 40:9;2", "SM 42:6;2", "SM 42:7;2", "SM 42:8;2", "SM 43:6;2", "SM 43:8;2", "SM 43:9;2",
              "SM 44:7;2", "SM 44:8;2", "SM 45:5;2|SM 15:2;2/30:3", "LPA 15:1", "LPA 16:1", "LPA 17:2", "LPA 20:5", "LPA 24:0", "LPC 23:0", "LPC O_18:4", "LPA 16:2", "LPA 17:0",
              "LPA 17:1", "LPA 22:3",
              "LPS 16:3", "PA 35:2|PA 17:1_18:1", "PC 29:2", "PC 34:6", "PC 52:0|PC 21:0_31:0", "PC O_35:5|PC O_18:5_17:0", "PC O_38:9", "PC O_42:6|PC O_26:6_16:0", "PC O_51:0",
               "PG 31:0|PG 14:0_17:0", "SM 41:10;2", "Cer 44:0;2|Cer 20:0;2/24:0", "Cer 44:3;4|Cer 28:2;3(FA 16:0)", "Cer 45:3;4|Cer 14:1;2/15:1;O(FA 16:0)", "Cer 45:4;4|Cer 27:2;3(FA 18:1)", "Cer 48:2;2|Cer 34:2;2/14:0",
                "Cer 48:3;2|Cer 18:1;2/30:2", "Cer 49:3;4|Cer 33:2;3(FA 16:0)", "Cer 50:3;2|Cer 18:1;2/32:2", "Cer 42:0;3|Cer 18:0;3/24:0", "DG 40:4|DG 18:1_22:3", "Hex2Cer 55:12;2", "HexCer 41:1;2|HexCer 18:1;2/23:0",
              "HexCer 55:14;2|HexCer 13:2;2/42:12", "HexCer 57:12;2|HexCer 13:2;2/44:10", "LPC 19:0_A", "LPC 19:0_B", "LPC 20:4", "LPC 22:1", "LPC 22:3", "LPC 28:0", "LPC O_18:4", "LPG 18:2",
              "LPG 20:3", "LPI 14:0", "LPS 20:3", "PC 24:0", "PC 27:1", "PC 28:4|PC 10:0_18:4", "PC 32:2", "PC 32:8", "PC 35:0|PC 17:0_18:0", "PC 36:3|PC 16:0_20:3", "PC 39:6|PC 19:2_20:4",
              "PC 40:10|PC 16:4_24:6", "PC 40:7|PC 16:1_24:6", "PC 42:6|PC 18:2_24:4", "PC 49:2|PC 17:0_32:2", "PC 52:11|PC 24:5_28:6", "PC O_34:8", "PC O_34:9", "PC O_36:7", "PC O_36:8", "PC O_36:9",
              "PC O_38:7", "PC O_38:8_A", "PC O_38:8_B", "PC O_38:9", "PC O_39:4", "PC O_39:8", "PC O_40:10", "PC O_40:11", "PC O_40:7", "PC O_40:8_A", "PC O_40:8_B", "PC O_40:9",
              "PC O_41:7", "PC O_41:8", "PC O_42:10", "PC O_42:7", "PC O_42:8", "PC O_42:9", "PC O_43:7", "PC O_44:10", "PC O_45:9", "PC O_47:0", "PC O_48:1", "PC O_48:7|PC O_24:6_24:1",
              "PC O_49:0", "PC O_51:0", "PC O_51:1", "PC O_53:1", "PE 36:10|PE 18:5_18:5", "PE 36:5", "PE 36:6", "PE O_19:0", "PE O_32:8", "PE O_36:5|PE O_18:4_18:1", "PE O_36:6|PE O_18:5_18:1",
              "PE O_40:1|PE O_22:0_18:1", "PE O_42:1|PE O_24:0_18:1", "PG 20:0|PG 10:0_10:0", "PI 35:0|PI 17:0_18:0", "PI 40:4|PI 18:1_22:3", "PI O_38:1|PI O_18:0_20:1", "PI O_50:7",
              "PS 17:0", "PS 32:4", "PS 32:6", "PS 34:0", "PS 34:5", "PS 34:6", "PS 36:3|PS 19:1_17:2", "PS 36:6", "PS 37:1", "PS 38:7|PS 20:3_18:4", "PS 40:1|PS 22:0_18:1", "PS 40:7",
                "SM 34:1;2_A", "SM 34:3;2|SM 18:1;2/16:2", "SM 40:1;2_A", "SM 40:2;2_A", "SM 40:2;2_B", "SM 40:2;2_C", "SM 41:0;2", "SM 41:10;2", "SM 50:0;2|SM 28:0;2/22:0",
              "SM 51:3;2", "SM 52:1;2|SM 33:0;2/19:1", "SM 56:2;2|SM 24:1;2/32:1", "SM 57:2;2", "SM 58:3;2|SM 18:1;2/40:2", "SM 59:3;2", "SM 60:10;2|SM 18:1;2/42:9", "SM 60:1;2|SM 18:1;2/42:0_A", "SM 60:1;2|SM 18:1;2/42:0_B",
              "SM 60:2;2|SM 18:1;2/42:1_A", "SM 60:2;2|SM 18:1;2/42:1_B", "SM 60:3;2|SM 18:1;2/42:2", "SM 61:3;2", "SM 62:1;2", "SM 65:0;2|SM 21:0;2/44:0", "SM 67:0;2|SM 34:0;2/33:0",
              "TG 67:2|TG 17:0_18:1_32:1", "TG 68:1|TG 16:0_18:0_34:1", "TG 72:2|TG 18:0_18:1_36:1", "TG 72:3|TG 18:1_18:1_36:1", "TG 72:3|TG 18:1_20:1_34:1", "TG 72:4|TG 18:1_18:1_36:2",
              "TG 74:2|TG 18:0_18:1_38:1", "TG 74:3|TG 18:1_18:1_38:1", "TG 74:4|TG 18:1_18:1_38:2", "CL 74:4|CL 18:0_18:1_18:0_20:3", "CL 72:5|CL 16:0_18:1_18:0_20:4", "CL 66:1|CL 16:0_16:0_16:0_18:1",
              "CL 66:3|CL 16:0_16:1_16:1_18:1_B", "CL 68:1|CL 16:0_18:0_16:0_18:1", "CL 70:1|CL 16:0_18:1_18:0_18:0", "CL 70:2|CL 16:0_18:0_18:1_18:1", "CL 70:2|CL 16:0_18:1_18:0_18:1",
              "CL 70:5|CL 16:0_18:1_16:0_20:4", "CL 72:3|CL 18:0_18:1_18:1_18:1", "CL 72:4|CL 18:0_18:2_18:1_18:1","CL 72:6|CL 16:1_20:3_18:1_18:1", "Cer 36:1;2|Cer 18:1;2/18:0_D", "Cer 59:2;4|Cer 41:1;3(FA 18:0)",
              "DG 35:2|DG 17:1_18:1","DG 46:2|DG 16:1_30:1", "FAHFA 24:2;O|FAHFA 16:1/8:0;O", "FAHFA 28:1;O|FAHFA 18:0/10:0;O", "FAHFA 34:3;O|FAHFA 16:1/18:1;O", "FAHFA 35:8;O|FAHFA 18:5/17:2;O",
              "FAHFA 36:3;O|FAHFA 18:1/18:1;O", "PC 32:3", "PC 34:3|PC 18:1_16:2", "PC 34:4", "PC 38:5", "PC 49:6", "PE 28:0|PE 14:0_14:0", "PI 38:4|PI 18:1_20:3", "PS 30:1|PS 14:0_16:1", "PC O_33:4|PC O_16:4_17:0",

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