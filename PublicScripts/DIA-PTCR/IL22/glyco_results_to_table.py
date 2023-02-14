import pandas as pd
import numpy as np
from copy import deepcopy
import os
import re


def translate_id(id):
    id = re.sub("\(Na\+\)", "", id)
    id = re.sub("\(Ca2\+\)", "", id)
    id = re.sub("\(Fe2\+\)", "", id)
    id = re.sub("\(2Na\+\)", "", id)
    try:
        masssearch = re.search("= (.*)m", id)
        mass = masssearch.group(1)
        mass = float(mass)
    except:
        mass = -1

    try:
        glycansearch = re.search("\+(.*)\)", id)
        glycan = glycansearch.group(1)
    except:
        glycan = None
    # print(id, mass, glycan)
    return mass, glycan


def match_to_glyco_df(gdf, glycan, mass):
    index = gdf.index[gdf["Glycan"] == glycan].to_numpy()
    if len(index) == 1:
        dbmass = float(gdf.iloc[index]["Monoisotopic mass"])
        # print(index, glycan, mass, dbmass, dbmass - mass)
        return index
    else:
        print("Did not match! Name:", glycan, "Mass?:", mass)
        return -1


os.chdir("C:\\Data\\Luis Genentech")

datafile = "SA15_LysC.xlsx"
datafile = "IL22_SA8_LysC.xlsx"
datafile = "IL22_SA4_tryp.xlsx"

glycanfile = "Complete Glycan List.csv"

outfile = datafile[:-5] + "_glycan_list.csv"
outfile2 = datafile[:-5] + "_glycan_list_abridged.csv"

df = pd.read_excel(datafile)

gdf = pd.read_csv(glycanfile)

print(df.keys())

df["Site"] = df["Site"].str.replace("~", "")

sites = np.unique(df["Site"])
print(sites)

foundsites = []
for n, s in enumerate(sites):
    if "N" in s:
        subdf = deepcopy(df[df["Site"] == s])
        sumint = np.sum(subdf["MS Area"])
        normint = np.array(subdf["MS Area"].to_numpy() / sumint * 100)
        subdf["Norm Area"] = normint
        gdf["S" + str(n + 1)] = 0
        for i, row in subdf.iterrows():
            mass, glycan = translate_id(row["Identification"])
            index = match_to_glyco_df(gdf, glycan, mass)
            if index >= 0:
                percent = row["Norm Area"]
                gdf.loc[index, "S" + str(n + 1)] += percent
        foundsites.append(s)

if len(foundsites) < 4:
    gdf["S4"] = 0
    index = gdf.index[gdf["Glycan"] == "Free"].to_numpy()[0]
    print(index)
    gdf.loc[index, "S4"] += 100

if len(foundsites) == 2:
    gdf["S3"] = 0
    index = gdf.index[gdf["Glycan"] == "Free"].to_numpy()[0]
    print(index)
    gdf.loc[index, "S3"] += 100

gdf.to_csv(outfile)

b1 = gdf["S1"] > 0
b2 = gdf["S2"] > 0
b3 = gdf["S3"] > 0
b4 = gdf["S4"] > 0

gdf2 = gdf[b1 + b2 + b3 + b4]

gdf2.to_csv(outfile2)
