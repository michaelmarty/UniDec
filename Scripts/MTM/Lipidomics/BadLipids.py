import os
import pandas as pd
from LipidFunctions import fix_name, fix_dash, fix_plasmalogens, msdialnames_to_LC

topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

badfile = r"Combined_Lipids_BadLipids.csv"
baddf = pd.read_csv(badfile)

# Import transition list
transfile = r"Combined_TransitionList.csv"
transdf = pd.read_csv(transfile)

transdf = msdialnames_to_LC(transdf)

# Cleanup names in transition list
for i, row in transdf.iterrows():
    name = row["Molecule Name"]
    lclass = row["Molecule List Name"]
    name = fix_dash(name)

    name, lclass = fix_name(name, lclass)
    if name.startswith("TAG ") and lclass == "TG":
        lclass = "TAG"
    if name.startswith("DAG ") and lclass == "DG":
        lclass = "DAG"
    if name.startswith("MAG ") and lclass == "MG":
        lclass = "MAG"

    name = fix_plasmalogens(name)
    transdf.at[i, "Molecule Name"] = name
    transdf.at[i, "Molecule List Name"] = lclass

# print(transdf["Molecule List Name"].unique())
# print(baddf["Molecule List Name"].unique())
#
# # Print unique adducts
# print(transdf["Precursor Adduct"].unique())
# print(baddf["Precursor Adduct"].unique())
#
# # Print all EtherMGDG in both with just the columns used for sorting
# print(transdf[transdf["Molecule List Name"] == "EtherMGDG"] [["Molecule List Name", "Molecule Name", "Precursor Adduct"]])
# print(baddf[baddf["Molecule List Name"] == "EtherMGDG"] [["Molecule List Name", "Molecule Name", "Precursor Adduct"]])

# Remove any in transdf that aren't in baddf
mergedf = pd.merge(transdf, baddf, on=["Molecule List Name", "Molecule Name", "Precursor Adduct"], how="inner")
# Sort by Molecule List Name and Molecule Name
mergedf = mergedf.sort_values(by=["Molecule List Name", "Molecule Name"])
print(mergedf["Molecule Name"].unique())

# Write to new file
badtransfile = transfile.replace(".csv", "_BadLipids.csv")
mergedf.to_csv(badtransfile, index=False)
