import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats
import pandas as pd
import os
import pandas as pd
import glob
import numpy as np

from openpyxl import workbook

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


def get_avg(mavg, mstd):
    mavg = np.mean(mutant)
    wavg = np.mean(wt)


R = 8.31446261815324 / 1000.
units = "kJ/mol"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"

temps = [15, 20, 25, 30, 35]

globaldir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\WT_W14A_"
fname = "Extract_grid_2D_Extract6.txt"

for t in temps:
    temp = t
    topdir = globaldir + str(t) + "C"
    print(topdir)
    files = []
    data = []

    for path in glob.glob(f'{topdir}/*/**/', recursive=True):
        testpath = os.path.join(path, fname)
        if os.path.isfile(testpath):
            files.append(testpath)
    print(files)
    exit()

    if len(files) < 1:
        print("Failed at t=", t)

    for f in files:
        data = np.loadtxt(f)
        data.append(alldata)
        print(alldata)


    wts = []
    mts = []
    maxl = 7

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

    print(mavg)

# np.savetxt("Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\2DCombiner.csv", data, delimiter=",", header="Lnum,X,Int")
