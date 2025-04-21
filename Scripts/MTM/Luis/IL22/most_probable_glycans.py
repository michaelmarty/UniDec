import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os



os.chdir("C:\Data\Luis Genentech\\MHC2")

file2 = "MHC2_DPA_tryp"
#file2 = "SA15_Trypsin"
file2 = file2 + "_glycan_list_abridged.csv"

df = pd.read_csv(file2)

print(df.keys())

dfs1 = df[df["S1"] > 0]
dfs2 = df[df["S2"] > 0]
dfs3 = df[df["S3"] > 0]
dfs4 = df[df["S4"] > 0]

lens = [len(dfs1), len(dfs2), len(dfs3), len(dfs4)]

indexes = np.array(list(np.ndindex(tuple(lens))))

probs = []
masses = []
s1s = []
s2s = []
s3s = []
s4s = []
for i in indexes:
    r1 = dfs1.iloc[i[0]].fillna(1)
    r2 = dfs2.iloc[i[1]].fillna(1)
    r3 = dfs3.iloc[i[2]].fillna(1)
    r4 = dfs4.iloc[i[3]].fillna(1)
    rs = [r1, r2, r3, r4]

    p = r1["S1"] / 100 * r2["S2"] / 100 * r3["S3"] / 100 * r4["S4"] / 100

    mass = np.sum([r["Monoisotopic mass"] for r in rs])

    probs.append(float(p))
    masses.append(mass)

    s1s.append(r1["Glycan"])
    s2s.append(r2["Glycan"])
    s3s.append(r3["Glycan"])
    s4s.append(r4["Glycan"])

newdf = pd.DataFrame()
newdf["Mass"] = masses
newdf["Probability"] = probs
newdf["Glycan1"] = s1s
newdf["Glycan2"] = s2s
newdf["Glycan3"] = s3s
newdf["Glycan4"] = s4s

newdf.to_excel(file2[:-4]+"_massnprob.xlsx")
