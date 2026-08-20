# Script to calculate DDG values from each quartet of extracted peak areas

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

#Definition of Kd and delta(deltaG)
def get_kd(mutant, wt, n=0, l=1, T=37):
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    TK = T + 273.15
    deltag = -R * TK * np.log(kd)
    return kd, deltag, n


R = 8.31446261815324 / 1000.
units = "kJ/mol"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"

# Data info
temps = [15, 20, 25, 30, 35]
globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final2\\WT_R224A Titrations All\\WT_R224A_"
fname = "Extract_grid_2D_Extract6.txt" # Sets the name of the file to be read in

#Selecting txt files
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
    if len(files)<2:
        print("Failed at t=", t)
        break
    kds = []
    wts = []
    mts = []
    maxl = 7

#Loading txt file data and normalization
    fig, ax = plt.subplots()
    for f in files:
        data = np.loadtxt(f)

        bm = data[:, 1] == 0
        mutant = data[bm][:, 2]
        wt = data[np.logical_not(bm)][:, 2]

        wt = wt[:maxl]
        mutant = mutant[:maxl]
        wt = wt / np.sum(wt)
        mutant = mutant / np.sum(mutant)
        wts.append(wt)
        mts.append(mutant)

        ax.plot(np.arange(len(wt)) + 0.5, wt, color="g")
        ax.plot(np.arange(len(mutant)), mutant, color="b")

        # Calculate the DDG for each
        kdlist = []
        for n in range(0, maxl - 1):
            for m in range(1, 2):# maxl - n):
                value = get_kd(mutant, wt, n, m, T=temp)
                kdlist.append(value)
        kds.append(kdlist)

    kds = np.array(kds)
    wts = np.array(wts)
    mts = np.array(mts)

    # Get Mean and Std. Dev
    wavg = np.mean(wts, axis=0)
    mavg = np.mean(mts, axis=0)

    wstd = np.std(wts, axis=0)
    mstd = np.std(mts, axis=0)

    #Printing Kd, ddG, ddG, n>n+1
    print(kds)

    # Make Plots
    maxn = maxl
    labels = ["0", "1", "2", "3", "4", "5", "6"]
    ax.bar(np.arange(maxn), mavg[:maxn], 0.45, color="b", yerr=mstd[:maxn], capsize=3, label="Mutant")
    ax.bar(np.arange(maxn) + 0.45, wavg[:maxn], 0.45, color="g", yerr=wstd[:maxn], capsize=3, label="WT")
    ax.set_xticks(np.arange(maxn) + 0.225, labels=labels[:maxn])
    ax.set_yticks([0, 0.5, 1.0], labels=["0", "%", "100"])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Number of Bound Lipids")
    ax.legend()
    ax.set_title(str(t))

    dgs = []

    for n in range(0, maxn - 1):
        avgdg = np.mean(kds[:, n, 1])
        sddg = np.std(kds[:, n, 1])

        dgvals = list(kds[:, n, 1])
        dgvals.append(avgdg)
        dgvals.append(sddg)

        dgs.append(dgvals)
        #ddG labels
        text = "\u0394\u0394G = " + str(round(avgdg, 2)) + " \u00B1 " + str(round(sddg, 2)) + " " + units
        ax.text(n, wavg[0] + 0.1 + 0.15 * n, text)
    # Save results
    np.savetxt(os.path.join(topdir, "DG_resultsN2_Fixed.txt"), dgs)
    plt.savefig(os.path.join(topdir, "PlotN2_Fixed.pdf"))
plt.show()
