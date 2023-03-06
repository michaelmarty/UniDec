import pandas as pd
import os

file = "C:\Data\Luis Genentech\Complete Glycan List PepFinder 1.0.csv"
file2 = "C:\Data\Luis Genentech\Complete Glycan List.csv"
#file = "C:\Data\Luis Genentech\human_glycans_formatted.csv"
#file2 = "C:\Data\Luis Genentech\human_glycans_formatted2.csv"

df = pd.read_csv(file)

nsas = []
ngals = []
ngns = []
nhexs = []
ants = ["Antenna 1", "Antenna 2", "Antenna 3", "Antenna 4"]
for index, row in df.iterrows():
    print(row)
    nsa = 0
    ngal = 0
    ngn = 0
    for a in ants:
        s = row[a]
        if "Gn-" in s:
            ngn += 1
        if "G-" in s:
            ngal += 1
        if "S-" in s:
            nsa += 1
    nman = row["Man"]-3
    nhex = nman + ngal

    nsas.append(nsa)
    ngals.append(ngal)
    ngns.append(ngn)
    nhexs.append(nhex)

df["NSa"] = nsas
df["NGn"] = ngns
df["NGal"] = ngals
df["NHex"] = nhexs

df.to_csv(file2)

print(df.keys())