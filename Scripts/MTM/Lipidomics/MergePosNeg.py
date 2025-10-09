import os
import pandas as pd
import numpy as np

def find_match_name(name, rt, refdf, rttol=0.1):
    newname = ""
    # Look for same name in negdf contained ahead of "|"
    rows = refdf[refdf["Metabolite name"].str.contains(name + "|", regex=False)]
    if len(rows) == 0:
        return newname
    if len(rows) == 1:
        refrt = rows.iloc[0]["Average Rt(min)"]
        refname = rows.iloc[0]["Metabolite name"]
        # If RTs are within 0.5 min, assume same lipid
        if abs(rt - refrt) < rttol:
            # print("Match found:", name, refname, rt, refrt)
            newname = refname
        # else:
        #     print("RTs not close:", name, rtpos, negname, negrt)
    else:
        # Find one with the closest rt
        rtdiff = np.abs(rows["Average Rt(min)"] - rt)
        rows = rows[rtdiff == rtdiff.min()]
        if len(rows) > 1:
            print("Warning Still multiple matches for", name, rt, rows[["Metabolite name", "Average Rt(min)"]])
        refrt = rows.iloc[0]["Average Rt(min)"]
        refname = rows.iloc[0]["Metabolite name"]
        if abs(rt - refrt) < 0.1:
            # print("Closest match found:", name, refname, rt, refrt)
            newname = refname
        # else:
        #     print("RTs not close:", name, rtpos, negname, negrt)
    return newname

def merge_pos_neg(posdf, negdf, rttol=0.1):
    print(posdf.columns)
    # Loop through posdf and look for names without "|"
    for i, row in posdf.iterrows():
        name = row["Metabolite name"]
        rtpos = row["Average Rt(min)"]
        if "|" not in name:
            newname = find_match_name(name, rtpos, negdf, rttol=rttol)
            if newname != "":
                posdf.at[i, "Metabolite name"] = newname
                print("Positive Updated name:", name, "to", newname)

    # Repeat in negative mode
    for i, row in negdf.iterrows():
        name = row["Metabolite name"]
        negrt = row["Average Rt(min)"]
        if "|" not in name:
            newname = find_match_name(name, negrt, posdf, rttol=rttol)
            if newname != "":
                negdf.at[i, "Metabolite name"] = newname
                print("Negative Updated name:", name, "to", newname)

    return posdf, negdf

if __name__ == "__main__":
    topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
    os.chdir(topdir)

    posfile = "PosIDsC1_Ref.csv"
    negfile = "NegIDsC1_Ref.csv"

    posdf = pd.read_csv(posfile)
    negdf = pd.read_csv(negfile)

    rttol = 0.1

    newposdf, newnegdf = merge_pos_neg(posdf, negdf, rttol)

    newposdf.to_csv("PosIDsC1_Ref.csv", index=False)
    newnegdf.to_csv("NegIDsC1_Ref.csv", index=False)

