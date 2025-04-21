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
    if glycan == "Unglycosylated":
        glycan = "Free"
    index = gdf.index[gdf["Glycan"] == glycan].to_numpy()
    if len(index) == 1:
        dbmass = float(gdf.iloc[index]["Monoisotopic mass"])
        # print(index, glycan, mass, dbmass, dbmass - mass)
        return index
    else:
        print("Did not match! Name:", glycan, "Mass?:", mass)
        return -1

def set_col_to_free(gdf, col):
    gdf[col] = 0
    index = gdf.index[gdf["Glycan"] == "Free"].to_numpy()[0]
    print("Setting column to all free:", col)
    gdf.loc[index, col] += 100


os.chdir("C:\\Data\\Luis Genentech\\Ovalbumin")

datafile = "ovalbumin_v8de.xlsx"
#datafile = "MHC2_DPA_v8de.xlsx"

glycanfile = "..\\Complete Glycan List.csv"

outfile = datafile[:-5] + "_glycan_list.csv"
outfile2 = datafile[:-5] + "_glycan_list_abridged.csv"

df = pd.read_excel(datafile)

gdf = pd.read_csv(glycanfile)

print(df.keys())

df["Site"] = df["Site"].str.replace("~", "")

sites = np.unique(df["Site"])
print(sites)

foundsites = []
snames = []
sindex = 1
for n, s in enumerate(sites):
    if "N" in s or "T" in s:
        subdf = deepcopy(df[df["Site"] == s])
        sumint = np.sum(subdf["MS Area"])
        normint = np.array(subdf["MS Area"].to_numpy() / sumint * 100)
        subdf["Norm Area"] = normint
        sname = "S" + str(sindex)
        gdf[sname] = 0
        for i, row in subdf.iterrows():
            mass, glycan = translate_id(row["Identification"])
            index = match_to_glyco_df(gdf, glycan, mass)
            if index >= 0:
                percent = row["Norm Area"]
                gdf.loc[index, sname] += percent
        foundsites.append(s)
        snames.append(sname)
        sindex += 1

if len(foundsites) < 4:
    set_col_to_free(gdf, "S4")

if len(foundsites) == 2:
    set_col_to_free(gdf, "S3")

print(snames)
bout =  gdf[snames[0]] > 0
for s in snames:
    b1 = gdf[s] > 0
    if not np.any(b1):
        set_col_to_free(gdf, s)
    bout = bout + b1


gdf.to_csv(outfile)

gdf2 = gdf[bout]

gdf2.to_csv(outfile2)
