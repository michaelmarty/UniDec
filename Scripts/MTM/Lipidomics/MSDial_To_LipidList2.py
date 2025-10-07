import pandas as pd
import os
from LipidFunctions import lipidlist_to_LC, check_FA_count

def clean_names(name):
    if "|" in name:
        name = name.split("|")[1]
    return name


def msdial_to_lipidlist(df):
    # Drop all where "Ontology" is Others or Unknown
    df = df[~df["Ontology"].isin(["Others", "Unknown"])]

    # Export Df with "Molecule List Name", "Molecule Name", "Precursor Adduct"
    df = df[["Ontology", "Metabolite name", "Adduct type"]]
    df = df.rename(columns={"Ontology": "Molecule List Name", "Metabolite name": "Molecule Name",
                            "Adduct type": "Precursor Adduct"})

    # Sort by Molecule List Name and Molecule Name
    df = df.sort_values(by=["Molecule List Name", "Molecule Name"])

    # Clean names
    df["Molecule Name"] = df["Molecule Name"].apply(clean_names)

    # Check FA count
    df = check_FA_count(df)

    # Convert to LC format
    df = lipidlist_to_LC(df)

    return df


topdir = "Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX\\FinalResults"

file = "PosIDsC1.csv"
os.chdir(topdir)

df = pd.read_csv(file)
df = msdial_to_lipidlist(df)
dfpos = df.copy()
print(df)
# Save to csv
df.to_csv("Pos_LipidList.csv", index=False)

file = "NegIDsC1.csv"
df = pd.read_csv(file)
df = msdial_to_lipidlist(df)
df.to_csv("Neg_LipidList.csv", index=False)
print(df)

# Merge pos and neg
df = pd.concat([dfpos, df])
df = df.drop_duplicates()

df = df.sort_values(by=["Molecule List Name", "Molecule Name"])
df.to_csv("Combined_Lipids.csv", index=False)
df.to_csv("Combined_LipidList.csv", columns=["Full Name Adduct"], index=False, header=False)
