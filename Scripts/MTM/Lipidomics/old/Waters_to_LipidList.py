import os
import pandas as pd
from LipidFunctions import lipidlist_to_LC, check_FA_count

def clean_names(name):
    # If _ + letter is at the end of the name, remove it
    if "_" in name and name.split("_")[-1].isalpha() and len(name.split("_")[-1]) == 1:
        name = "_".join(name.split("_")[:-1])
    # If / + letter is at the end of the name, remove it
    return name

def waters_to_lipidlist(df):
    # Keep only columns "Molecule List Name" and "Name" and "Precursor Adduct"
    df = df[["Molecule List Name", "Name", "Precursor Adduct"]]

    # Rename "Name" to "Molecule Name"
    df = df.rename(columns={"Name": "Molecule Name"})

    # Clean names
    df["Molecule Name"] = df["Molecule Name"].apply(clean_names)

    # Sort by Molecule List Name and Molecule Name
    df = df.sort_values(by=["Molecule List Name", "Molecule Name"])

    # Check FA count
    df = check_FA_count(df)

    # Convert to LC format
    df = lipidlist_to_LC(df)
    return df

topdir = r"Z:\Group Share\Annika\Stellar"
os.chdir(topdir)
file = "Waters_Brain_Both_output.xlsx"

df = pd.read_excel(os.path.join(topdir, file))
df = waters_to_lipidlist(df)
print(df)
dfbrain = df.copy()

file2 = "Waters_Ecoli_output.xlsx"
df = pd.read_excel(os.path.join(topdir, file2))
# Add column of "Precursor Adduct" with value "[M-H]"
df["Precursor Adduct"] = "[M-H]"
df = waters_to_lipidlist(df)
print(df)


# Merge pos and neg
df = pd.concat([dfbrain, df])
df = df.drop_duplicates()
df = df.sort_values(by=["Molecule List Name", "Molecule Name"])
df.to_csv("Waters_Combined_LipidList.csv", index=False)