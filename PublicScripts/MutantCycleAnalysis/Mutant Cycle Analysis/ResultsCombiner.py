import pandas as pd
import os

"""
path = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final4\\"
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

topdir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\"
os.chdir(topdir)
temps = [15, 20, 25, 30, 35]

dfs = []

edges = [("WT", "R224A"), ("WT", "W14A"), ("WT", "DM2"), ("W14A", "DM2"), ("R224A", "DM2"), ("R224A", "W14A")]
pairs = ["WT_R224A", "WT_W14A", "WT_DM2", "W14A_DM2", "R224A_DM2", "R224A_W14A"]
for ed in range(0, 6):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2+" Titrations All"
    pair = pairs[ed]
    path = os.path.join(topdir, dirname)

    file = "AnalysisResults2.csv"
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
df.to_excel("CombinedResults2.xlsx")
print(df)