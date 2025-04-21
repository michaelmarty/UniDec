import pandas as pd
import os

'''
1: Molecule List Name
2: Moleule Name
3: Molecule Formula
4: Precursor Adduct
5: Precursor Charge
6: Frament Ion Names
7: Product Mz
8: Product CHarge
9: Product Ion Formula
10: Product Adduct
'''

toppath = "Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\All_Reps"
os.chdir(toppath)

file = "Height_1_2024_08_23_14_32_45.txt"

df = pd.read_table(file, header=4)

df = df[df["MS/MS matched"]]

simp_names = []
full_names = []
for i, row in df.iterrows():
    mol_list_name = row["Ontology"]
    mol_name = row["Metabolite name"]
    if "|" in mol_name:
        names = mol_name.split("|")
        simp_name = names[0]
        full_name = names[1]
    else:
        simp_name, full_name = mol_name, mol_name
    simp_names.append(simp_name)
    full_names.append(full_name)

    mol_formula = row["Formula"]
    prec_adduct = row["Adduct type"].split("]")[0] + "]"
    prec_mz = row["Average Mz"]
    prec_charge = row["Adduct type"].split("]")[1]

df["Simp Name"] = simp_names
df["Full Name"] = full_names
df = df[df["Full Name"].str.contains("RIKEN") == False]

df = df.sort_values("Full Name")
print(df["Full Name"])
df.to_csv("Filtered.csv")
df.to_csv("LipidList.csv", columns=["Full Name"], index=False, header=False)
