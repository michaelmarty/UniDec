import numpy as np
import pandas as pd
import scipy.stats as stats
from copy import deepcopy
import os

def correct_data(df):
    print(df)

    n = 16
    dof = n - 1
    t = stats.t.ppf(0.975, dof)

    # Correct Occupancy
    std = df["Std Occupancy"].to_numpy()
    ci = t * std / np.sqrt(n)
    df["CI"] = ci

    lower_bound = df["Average Occupancy"] - df["CI"]
    b1 = lower_bound < 0

    occ = df["Average Occupancy"].to_numpy()
    occ = deepcopy(occ)
    occ[b1] = 0

    df["Corrected Occupancy"] = occ

    # Correct Residence Time
    rstd = df["Std Residence Time"].to_numpy()
    ci = t * rstd / np.sqrt(n)
    df["Res CI"] = ci

    lower_bound_res = df["Average Residence Time"] - df["Res CI"]
    b2 = lower_bound_res < 0

    res = df["Average Residence Time"].to_numpy()
    res = deepcopy(res)

    res[b2] = 0

    df["Corrected Residence Time"] = res
    return df

dir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis\\VMD\\CG_Figures\\Mutants_Dataset\\W14A"
os.chdir(dir)

lipids = ["CDL2", "POPE", "POPG"]

for l in lipids:
    file = "Averaged_Dataset_" + l + "_W14A.xlsx"
    outfile = "Averaged_Dataset_" + l + "_W14A_Corrected.xlsx"

    df = pd.read_excel(file, index_col=0)

    df = correct_data(df)

    df.to_excel(outfile)

