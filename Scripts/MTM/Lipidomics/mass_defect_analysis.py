import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

kendrick_ref_mass = 14.01565
mass_h = 1.007276467

def kendrick_analysis(mz, charge):
    mass = np.abs(mz * charge) - charge * mass_h
    kmass = mass * 14 / kendrick_ref_mass
    kmd = kmass - np.round(kmass)
    return mass, kmass, kmd

"""
toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

file = "AreasAllFixed.csv"
fulldf = pd.read_csv(file)
fulldf = fulldf[fulldf["MS/MS matched"]==False]
fullmz = fulldf["Average Mz"].to_numpy()
fullz = fulldf.Charge.to_numpy()

df = pd.read_csv("FilteredFixedRT.csv")
print(df.keys())

b1 = df.Ontology == "PG"
b2 = df.Ontology == "PE"
b3 = df.Ontology == "CL"
df = df[b3]

refmz = df["Reference m/z"].to_numpy()
mz = df["Average Mz"].to_numpy()
rt = df["Average Rt(min)"].to_numpy()
ontcodes = df["Ontology"].astype("category").cat.codes
charge = df.Charge.to_numpy()
print(np.unique(charge))

mass, kmass, kmd = kendrick_analysis(mz, charge)
refmass, refkmass, refkmd = kendrick_analysis(refmz, charge)
fmass, fkmass, fmd = kendrick_analysis(fullmz, fullz)
deltamd = np.abs(kmd - refkmd)

b2 = deltamd > 0.005

print(df[b2]["Full Name"])

plt.scatter(mass, kmd, marker="o", c=ontcodes)
plt.scatter(refmass, refkmd, marker="o", c="green")
# plt.scatter(kmass[b2], kmd[b2], marker="o", c="red")

buffer = 0.1
kmassrange = [np.amin(kmass) * (1 - buffer), np.amax(kmass) * (1 + buffer)]
kmdrange = [np.amin(kmd) + np.amin(kmd) * buffer, np.amax(kmd) + np.amax(kmd) * buffer]


b1 = (fkmass < kmassrange[1]) * (fkmass > kmassrange[0]) * (fmd < kmdrange[1]) * (fmd > kmdrange[0])



#plt.scatter(fmass[b1], fmd[b1], c="red")

plt.show()
"""

