import pandas as pd
import os
import csv

"""
path = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final9\\"
os.chdir(path)
temps = [15, 20, 25, 30, 35]
pairs = ["WT_R224A", "WT_W14A", "WT_DM2", "W14A_DM2", "R224A_DM2"]

dfs = []
for i in range(0, len(pairs)):
    pair = pairs[i]
    file = pair+".csv"
    df = pd.read_csv(os.path.join(path, file))
    pboth = [pair for n in range(len(df))]
    p1, p2 = pair.split("_")
    p1 = [p1 for n in range(len(df))]
    p2 = [p2 for n in range(len(df))]
    df["Mutant"]=pboth
    df["P1"] = p1
    df["P2"] = p2
    dfs.append(df)

df = pd.concat(dfs)
df.to_excel("CombinedResults.xlsx")
print(df)"""

topdir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis2\\POPE\\"
os.chdir(topdir)
temps = [15, 20, 25, 30, 35]

dfs = []

edges = [("WT", "W14A"), ("WT", "W206A"), ("WT", "W209A"), ("WT", "F10A"), ("WT", "F13A"), ("WT", "F196A"), ("WT", "Y100A"), ("WT", "R3A"), ("WT", "R75A"), ("WT", "R224A"), ("WT", "K4A"), ("WT", "K79A"), ("WT", "K155A")]
pairs = ["WT_W14A", "WT_W206A", "WT_W209A", "WT_F10A", "WT_F13A", "WT_F196A", "WT_Y100A", "WT_R3A", "WT_R75A", "WT_R224A", "WT_K4A", "WT_K79A", "WT_K155A"]
for ed in range(0, 13):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    dirname = e1 + "_" + e2+" Titrations All"
    pair = pairs[ed]
    path = os.path.join(topdir, dirname)

    file = "AnalysisResults2.csv"
    df = pd.read_csv(os.path.join(path, file), nrows=3)
    pboth = [pair for n in range(len(df))]
    p1, p2 = pair.split("_")
    p1 = [p1 for n in range(len(df))]
    p2 = [p2 for n in range(len(df))]
    df["Mutant"]=pboth
    df["P1"] = p1
    df["P2"] = p2
    dfs.append(df)

df = pd.concat(dfs)
df.to_excel("CombinedResults_POPE_test.xlsx")
print(df)