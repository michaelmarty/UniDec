import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import pandas as pd


temps = [15, 20, 25, 30, 35]
globaldir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\WT_W14A_"
fname = "Extract_grid_2D_Extract6.txt"


dflist = []


for t in temps:
    temp = t
    topdir = globaldir+str(t)+"C"
    print(topdir)
    files = []
    for path in glob.glob(f'{topdir}/*/**/', recursive=True):
        testpath = os.path.join(path,fname)
        if os.path.isfile(testpath):
            files.append(testpath)
    print(files)
    if len(files)<1:
        print("Failed at t=", t)
        break
    kds = []
    wts = []
    mts = []
    maxl = 7


    for f in files:
        data = np.loadtxt(f)

        bm = data[:, 1] == 0
        mutant = data[bm][:, 2]
        wt = data[np.logical_not(bm)][:, 2]

        wt = wt[:maxl]
        mutant = mutant[:maxl]

        wts.append(wt)
        mts.append(mutant)

    wts = np.array(wts)
    mts = np.array(mts)

    wavg = np.mean(wts, axis=0)
    mavg = np.mean(mts, axis=0)

    wstd = np.std(wts, axis=0)
    mstd = np.std(mts, axis=0)

    df = pd.DataFrame()

    df['Mutant Average'] = pd.Series(mavg)
    df['Mutant SD'] = pd.Series(mstd)
    df['WT Average'] = pd.Series(wavg)
    df['WT SD'] = pd.Series(wstd)

    dflist.append(df)

output_df = pd.concat(dflist)
print(output_df)

output_df.to_csv("Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\2DCombiner.csv")
