import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 1


# rcParams['font.size']=8

def get_kd(mutant, wt, n=0, l=1, T=37):
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    TK = T + 273.15
    deltag = -R * TK * np.log(kd)
    return kd, deltag, deltag / l, n, l, n + l


R = 8.31446261815324 / 1000.
units = "kJ/mol"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"

temps = [15, 20, 25, 30, 35]
globaldir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_R224A Titrations All\\WT_R224A_"
fname = "Extract_grid_2D_Extract6.txt"

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
        #print("data")
        #print(data)

        bm = data[:, 1] == 0
        print(bm)
        mutant = data[bm][:, 2]
        wt = data[np.logical_not(bm)][:, 2]

        wt = wt[:maxl]
        mutant = mutant[:maxl]

        wts.append(wt)
        mts.append(mutant)
        #print(mts)

        #value = get_kd(mutant, wt, 0, maxl - 1, T=temp)
        #print("Overall", value)
        kdlist = []

        for n in range(0, maxl - 1):
            for m in range(1, 2):# maxl - n):
                value = get_kd(mutant, wt, n, m, T=temp)
                kdlist.append(value)
        kds.append(kdlist)

    kds = np.array(kds)
    wts = np.array(wts)
    mts = np.array(mts)

    wavg = np.mean(wts, axis=0)
    mavg = np.mean(mts, axis=0)

    wstd = np.std(wts, axis=0)
    mstd = np.std(mts, axis=0)

    print(mavg)
    print(mstd)
    print(wavg)
    print(wstd)

#np.savetxt("Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\2DCombiner.csv", mavg, mstd, wavg, wstd, delimiter=",", header="mavg, mstd, wavg, wstd")

