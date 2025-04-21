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

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

file = "Area_0_20227291548.txt"
file = "MsmsIncluded_0_2022729127.txt"
file = "RT_0_2022861815.txt"

df = pd.read_table(file, header=4)

df = df[df["MS/MS matched"]]

simp_names = []
full_names = []
full_names_adduct = []
charges = []
for i, row in df.iterrows():
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

df["Simp Name"] = simp_names
df["Full Name"] = full_names
df["Full Name Adduct"] = full_names_adduct
df["Charge"] = charges
df = df[df["Full Name"].str.contains("RIKEN") == False]
df = df.sort_values("Full Name")
df.to_csv("FilteredFixedRT.csv")
df.to_csv("LipidList.csv", columns=["Full Name Adduct"], index=False, header=False)
