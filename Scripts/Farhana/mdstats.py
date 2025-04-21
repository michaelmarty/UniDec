import pandas as pd
import numpy as np
import os
from scipy import stats

file = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis\\VMD\\CG_Figures\\MutantsResTimeAll.xlsx"

os.chdir(os.path.dirname(file))

outfile = "MDStats.xlsx"

muts = ["WT", "R3A", "W14A", "F196A", "W206A"]

comparisons = [["WT", "R3A"], ["WT", "W14A"], ["WT", "F196A"], ["WT", "W206A"]]

dfout = pd.DataFrame()

n = 16
dof = n - 1
t = stats.t.ppf(0.975, dof)
print(t, t / np.sqrt(n))

for c in comparisons:
    df1 = pd.read_excel(file, sheet_name=c[0], index_col=0)
    df2 = pd.read_excel(file, sheet_name=c[1], index_col=0)

    meanname = "Average Residence Time"
    stdname = "Std Residence Time"

    diff = df1[meanname] - df2[meanname]
    diff_std = np.sqrt(df1[stdname] ** 2 + df2[stdname] ** 2)

    ci = t * diff_std / np.sqrt(n)

    # check if ci crosses 0, if so set diff to 0

    diff = diff.to_numpy()
    diff = diff.copy()

    bpos = diff > 0
    bneg = diff < 0

    lower_bound = diff - ci
    b1 = lower_bound < 0
    b1 = b1 & bpos
    diff[b1] = 0

    upper_bound = diff + ci
    b2 = upper_bound > 0
    b2 = b2 & bneg
    diff[b2] = 0

    dfout["Residue"] = df1["Residue"]

    dfout[c[0] + " - " + c[1]] = diff
    dfout[c[0] + " - " + c[1] + " CI"] = ci
    dfout[c[0] + " - " + c[1] + " STD"] = diff_std

print(dfout)

dfout.to_excel(outfile)