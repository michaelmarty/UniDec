import pandas as pd
import numpy as np
import os
import re


def get_lipid_array(name):
    sorts = np.sort(re.split('_| ', name))
    out = ""
    for s in sorts:
        out = out + s
    return out


def MSMS_to_array(string):
    a1 = string.split(" ")
    spec = [a.split(":") for a in a1]
    return np.array(spec).astype(float)


adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"]])

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

df = pd.read_csv("Filtered.csv")
sorts = []
for i, row in df.iterrows():
    try:
        sorts.append(get_lipid_array(row["Full Name"]))
    except Exception as e:
        print(e)
        print(row["Full Name"])
        sorts.append("")
df["Sorts"] = sorts

tl = pd.read_excel("LC_transitions.xls")
print(tl.keys())

good = []
rts = []
for i, row in tl.iterrows():

    # Match Name
    rname = row["PrecursorName"]
    rname2 = rname.replace("-", "_")
    try:
        sname = get_lipid_array(rname2)
    except Exception as e:
        print(e)
        print(rname2)
        sname = "badbadbad"
    subdf = df[df["Sorts"] == sname]

    # Match Charge
    if len(subdf) > 0:
        charge = row["PrecursorCharge"]
        subdf = subdf[subdf["Charge"] == charge]

    # Match to Spectrum
    if len(subdf) > 0:
        match = subdf.iloc[0]
        rts.append(match["Average Rt(min)"])
        if row["ProductName"] == "precursor":
            good.append(True)
        else:
            testmz = row["ProductMz"]
            spectrum = MSMS_to_array(str(match["MS/MS spectrum"]))
            mindiff = np.amin(np.abs(testmz - spectrum[:, 0]))
            ppmerror = mindiff / testmz * 1e6

            if ppmerror < 50:
                good.append(True)
            else:
                good.append(False)

    else:
        good.append(False)
        rts.append(0)
        if row["ProductName"] == "precursor":
            print(rname)
            pass

tl["Explicit Retention Time"] = rts
tl["Explicit Retention Time Window"] = [0.5 for i in range(len(tl))]
tl2 = tl[good]
print(len(tl))
print(len(tl2))
tl2.to_csv("Filtered_Transitions.csv")
