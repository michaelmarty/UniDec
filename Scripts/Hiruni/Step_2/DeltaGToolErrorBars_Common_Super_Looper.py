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

temps = [15, 20, 25, 30, 35, 40, 45]
fname = "Extract_grid_2D_Extract1.txt"
globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final\\"

os.chdir(globaldir)
edges = [("WT", "R224A"), ("WT", "R75A"), ("WT", "DM"), ("R224A", "DM"), ("R75A", "DM")]
for ed in range(0, 5):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    pathtop = os.path.join(globaldir, dirname)
    print(os.path.isdir(pathtop), pathtop)
    out = True

    temps = [15, 20, 25, 30, 35, 40]
    maxn = 6
    alldata = []
    avgstd = []
    reps = []
    for t in temps:
        newpath = os.path.join(pathtop, e1 + "_" + e2 + "_" + str(t) + "C")
        os.chdir(newpath)
        #data = np.loadtxt("DG_resultsN2_Fixed.txt")
        temp = t
        topdir = newpath
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

        #print(kds)

        maxn = maxl
        labels = ["0", "1", "2", "3", "4", "5", "6"]
        ax.bar(np.arange(maxn), mavg[:maxn], 0.45, color="b", yerr=mstd[:maxn], capsize=3, label="Mutant")
        ax.bar(np.arange(maxn) + 0.45, wavg[:maxn], 0.45, color="g", yerr=wstd[:maxn], capsize=3, label="WT")
        ax.set_xticks(np.arange(maxn) + 0.225, labels=labels[:maxn])
        ax.set_yticks([0, 0.5, 1], labels=["0", "%", "100"])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Number of Bound Lipids")
        ax.legend()
        ax.set_title(str(t))

        dgs = []

        #print('AVG DG:', np.mean(kds[:, :, 2], axis=0))
        #print('STD DG:', np.std(kds[:, :, 2], axis=0))
        for n in range(0, maxn - 1):
            #b1 = kds[0, :, 3] == n
            b1 = n
            avgdg = np.mean(kds[:, b1, 2])
            sddg = np.std(kds[:, b1, 2])

            dgvals = list(kds[:, n, 2])
            dgvals.append(avgdg)
            dgvals.append(sddg)

            dgs.append(dgvals)

            text = "\u0394\u0394G = " + str(round(avgdg, 2)) + " \u00B1 " + str(round(sddg, 2)) + " " + units
            ax.text(n, wavg[0] + 0.2 + 0.05 * n, text)

        np.savetxt(os.path.join(topdir, "DG_resultsN2_Fixed_2.txt"), dgs)
        plt.savefig(os.path.join(topdir, "PlotN2_Fixed.pdf"))
plt.show()
