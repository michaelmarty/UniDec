import pandas as pd
import os
from LipidFunctions import lipidlist_to_LC

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

    # Convert to LC format
    df = lipidlist_to_LC(df)

    return df

topdir = r"Z:\Group Share\Annika\Stellar"

file = "Sciex_Pos.csv"
os.chdir(topdir)

df = pd.read_csv(file)
df = msdial_to_lipidlist(df)
dfpos = df.copy()
print(df)
# Save to csv
df.to_csv("Sciex_Pos_LipidList.csv", index=False)

file = "Sciex_Neg.csv"
df = pd.read_csv(file)
df = msdial_to_lipidlist(df)
df.to_csv("Sciex_Neg_LipidList.csv", index=False)
print(df)

# Merge pos and neg
df = pd.concat([dfpos, df])
df = df.drop_duplicates()
df = df.sort_values(by=["Molecule List Name", "Molecule Name"])
df.to_csv("Sciex_Combined_LipidList.csv", index=False)
