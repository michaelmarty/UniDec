import pandas as pd
import os
from LipidFunctions import msdial_to_lipidlist

# First run AddRefMSMS.py to add reference spectra to the msdial files
# Then run this to convert to LipidList format
# Then run MSDial_to_TransitionList.py to make transition list
# Then Run LCexe.py
# Then Run BadLipids.py to filter out bad lipids


topdir = "Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX\\FinalResults"

file = "PosIDsC1_Ref.csv"
os.chdir(topdir)

df = pd.read_csv(file)
df = msdial_to_lipidlist(df)
dfpos = df.copy()
print(df)
# Save to csv
df.to_csv("Pos_LipidList.csv", index=False)

file = "NegIDsC1_Ref.csv"
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
