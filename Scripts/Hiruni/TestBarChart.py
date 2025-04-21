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


def Kfunc(t, h, s):
    return np.exp((-(h / (R * t)) + s / R))


topdir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\WT_K4A\\"

os.chdir(topdir)
edges = [("WT", "K4A")]
for ed in range(0, 1):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    path = os.path.join(topdir, dirname)
    print(os.path.isdir(path), path)
    out = True

    temps = [25, 37]
    maxn = 6
    alldata = []
    avgstd = []
    for t in temps:
        newpath = os.path.join(path, e1 + "_" + e2 + "_" + str(t) + "C")
        os.chdir(newpath)
        data = np.loadtxt("DG_resultsN2_Fixed.txt")

        for i, d in enumerate(data):
            avgstd.append([i, t, d[-2], d[-1], 0, 0, 0, 0, 0, 0, 0, 0])
            for j, x in enumerate(d[:-2]):
                alldata.append([i, t, x, j])
    os.chdir(path)
    alldata = np.array(alldata)
    avgstd = np.array(avgstd)

        plt.figure()
                x = []
                y = []
                e = []
                for l in range(maxn):
                    b1 = avgstd[:, 0] == l
                    meanstd = avgstd[b1][2]
                    x.append(l + 1)
                    y.append(meanstd[2])  # 10
                    e.append(meanstd[3])  # 11

                plt.bar(x, y, yerr=e)  # , tick_labels=x+1)
                plt.gca().set_ylim(-0.5, 1.0)
                plt.yticks([-0.5, 0, 0.5, 1])
                plt.xlabel("Number of Lipids")
                plt.ylabel("Delta G")
                if out:
                    plt.savefig(os.path.join(topdir, "VH_Plot2.pdf"))
                # plt.show()