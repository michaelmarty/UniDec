import pandas as pd
import numpy as np
import os

topdir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis\\VMD\\CG_Figures\\Mutants_Dataset\\R3A"

os.chdir(topdir)


def average_df(df):
    avg_occ = []
    std_occ = []
    avg_restime = []
    std_restime = []
    uresidues = df["Residue"].unique()

    for res in uresidues:
        res_df = df[df["Residue"] == res]
        occ = res_df["Occupancy"]
        restime = res_df["Residence Time"]

        avg_occ.append(occ.mean())
        std_occ.append(occ.std())
        avg_restime.append(restime.mean())
        std_restime.append(restime.std())

    newdf = pd.DataFrame({"Residue": uresidues, "Average Occupancy": avg_occ, "Std Occupancy": std_occ,
                          "Average Residence Time": avg_restime, "Std Residence Time": std_restime})
    return newdf

def merge_replicates(df):
    uresidues = df["Residue"].unique()

    avg_occ = []
    std_occ = []
    avg_restime = []
    std_restime = []

    for res in uresidues:
        res_df = df[df["Residue"] == res]
        occ = res_df["Average Occupancy"]
        restime = res_df["Average Residence Time"]

        ostd = res_df["Std Occupancy"]
        rstd = res_df["Std Residence Time"]

        avg_occ.append(occ.mean())
        std_occ.append(np.sqrt(ostd.pow(2).sum())/len(ostd))
        avg_restime.append(restime.mean())
        std_restime.append(np.sqrt(rstd.pow(2).sum())/len(rstd))

    newdf = pd.DataFrame({"Residue": uresidues, "Average Occupancy": avg_occ, "Std Occupancy": std_occ,
                            "Average Residence Time": avg_restime, "Std Residence Time": std_restime})

    return newdf


edges = [("CDL2", "R3A"), ("POPE", "R3A"), ("POPG", "R3A")]

for ed in range(len(edges)):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    print(e1, e2)
    dfs = []
    rep = [1, 2, 3, 4]
    for r in rep:
        df = pd.read_csv("Dataset_" + e1 + "_" + e2 + "_" + str(r) + ".csv")
        avg_df = average_df(df)
        dfs.append(avg_df)
        print(dfs)

    avg_df = merge_replicates(pd.concat(dfs))
    avg_df.to_excel("Averaged_Dataset_" + e1 + "_" + e2 + ".xlsx")

