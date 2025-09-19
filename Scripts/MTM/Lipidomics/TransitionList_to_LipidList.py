import pandas as pd
import os

from LipidFunctions import lipidlist_to_LC

topdir = r"Z:\Group Share\Annika\Stellar"

file = "Thermo Merged Transition List.csv"
os.chdir(topdir)

df = pd.read_csv(file)

print(len(df))
# Filter for only unique Molecule Name
df = df.drop_duplicates(subset=["Molecule Name"])
print(len(df["Molecule Name"]))

# Export Df with "Molecule List Name", "Molecule Name", "Precursor Adduct"
df = df[["Molecule List Name", "Molecule Name", "Precursor Adduct"]]

df = lipidlist_to_LC(df)
# Save to csv
df.to_csv("Thermo Merged LipidList.csv", index=False)