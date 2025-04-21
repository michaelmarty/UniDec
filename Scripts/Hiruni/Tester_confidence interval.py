import pandas as pd
import os

"""
path = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final\\Analysis"
os.chdir(path)
temps = [15, 20, 25, 30, 35, 40]
pairs = ["WT_R224A", "WT_R75A", "WT_DM", "R75A_DM", "R224A_DM"]

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

