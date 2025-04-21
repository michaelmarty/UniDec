import pandas as pd
import numpy as np
import os

adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"]])

toppath = "D:\\Johnny\\Lipidomics\\Centroid"
os.chdir(toppath)

file = "Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\All_Reps\\Height_1_2024_08_23_14_32_45.txt"

df = pd.read_table(file)
df.columns = df.iloc[3]
df = df.iloc[4:]

to_remove = []

simp_names = []
full_names = []
full_names_adduct = []
charges = []
for i, row in df.iterrows():
    if row["m/z matched"] == False:
        to_remove.append(i)
    mol_list_name = row["Ontology"]
    mol_name = row["Metabolite name"]
    prec_adduct = row["Adduct type"]
    print(prec_adduct)
    if prec_adduct in adduct_translator[:, 0]:
        index = np.where(adduct_translator[:, 0] == prec_adduct)[0][0]
        adduct_row = adduct_translator[index]
        aname = adduct_row[-1]
        charge = adduct_row[-2]
    else:
        aname = ""
        charge = 0
    if "|" in mol_name:
        names = mol_name.split("|")
        simp_name = names[0]
        full_name = names[1]
    else:
        simp_name, full_name = mol_name, mol_name
    simp_names.append(simp_name)
    full_names.append(full_name)
    charges.append(charge)
    full_names_adduct.append(str(full_name) + str(aname))

to_remove_arr = np.array(to_remove)

np.sort(to_remove_arr)[::-1]



df["Simp Name"] = simp_names
df["Full Name"] = full_names
df["Full Name Adduct"] = full_names_adduct
df["Charge"] = charges

for i in to_remove_arr:
    df.drop(i)

df = df[df["Full Name"].str.contains("RIKEN") == False]
df = df.sort_values("Full Name")
df.to_csv("FilteredFixedRT.csv")
df.to_csv("LipidList.csv", columns=["Full Name Adduct"], index=False, header=False)
