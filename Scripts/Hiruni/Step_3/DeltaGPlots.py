import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.75
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 16
rcParams['font.sans-serif'] = "Arial"

R = 8.31446261815324 / 1000.
units = "kJ/mol"

#Selecting the file
topdir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final2\\"

os.chdir(topdir)
edges = [("WT", "R224A")]
for ed in range(0, 1):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    # if e2 is "DM":
    # e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    path = os.path.join(topdir, dirname)
    print(os.path.isdir(path), path)
    out = True

    #Setting the Temperature
    t = 35
    maxn = 6
    temp_val = 35

    alldata = []
    avgstd = []

    newpath = os.path.join(path, e1 + "_" + e2 + "_" + str(t) + "C")
    os.chdir(newpath)
    data = np.loadtxt("DG_resultsN2_Fixed.txt")

    for i, d in enumerate(data):
        avgstd.append([i, t, d[-2], d[-1], 0, 0])
        for j, x in enumerate(d[:-2]):
            alldata.append([i, t, x, j])


    os.chdir(path)
    alldata = np.array(alldata)
    avgstd = np.array(avgstd)

    for l in range(maxn):
        lipidn = l
        b1 = alldata[:, 0] == lipidn
        data = alldata[b1]
        b1 = avgstd[:, 0] == l
        meanstd = avgstd[b1]

        df = len(data) - 1
        tval = abs(stats.t.ppf(0.05 / 2, df))

        dGAvg = np.mean(data[:,2])
        dGStd = (np.std(data[:,2])) * (tval / np.sqrt(len(data)))

        b1 = avgstd[:, 0] == lipidn
        avgstd[b1, -2] = dGAvg
        avgstd[b1, -1] = dGStd
    print(avgstd)

    if out:
        np.savetxt("AnalysisResults3.csv", avgstd, delimiter=",", header="Lnum,Temp,G,GE, GA,GAE")

    # plt.tight_layout()
    # plt.show()
    if True:

        plt.figure()
        x = []
        y = []
        e = []
        for l in range(maxn):
            b1 = avgstd[:, 0] == l
            meanstd = avgstd[b1][0]
            x.append(l + 1)
            y.append(meanstd[-2])  # 10 # 2 means avg ddG from data and 10 means avg ddG from fit
            e.append(meanstd[-1])  # 11 # or 3 for std ddG from data and 11 for std ddG from fit

        plt.bar(x, y, yerr=np.abs(e))  # , tick_labels=x+1)
        plt.gca().set_ylim(-0.5, 1.0)
        plt.yticks([-0.5, 0, 0.5, 1.0])
        plt.xlabel("Number of Lipids")
        plt.ylabel("Delta G")
        if out:
            plt.savefig(os.path.join(topdir, "VH_Plot2.pdf"))
        # plt.show()

plt.show()
