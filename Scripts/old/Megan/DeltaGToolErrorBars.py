import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 1


# rcParams['font.size']=8

def get_kd(mutant, wt, n=0, l=1):
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    deltag = -R * T * np.log(kd)
    return kd, deltag, deltag / l


R = 8.31446261815324 / 1000.
units = "kJ/mol"
T = 303.15
#R = 1.98720425864083 / 1000.
#units = "kcal/mol"

topdir = "Z:\\Group Share\\Megan E\\AqpZ\\kD_analysis\\R75A_R224A_titrations\\all_CL_reps\\"

files = ["rep1\\AqpZ_R75A_R224A_80CL_unidecfiles", "rep2\\AqpZ_R75A_R224A_80CL_unidecfiles", "rep3\\AqpZ_R75A_R224A_80CL_unidecfiles"]

#files = ["rep1\\AqpZ_R75A_R224A_160CL_unidecfiles", "rep2\\AqpZ_R75A_R224A_160CL_unidecfiles"]#, "rep3\\AqpZ_R75A_R224A_160CL_unidecfiles"]
'''
topdir = "Z:\\Group Share\\Megan E\\AqpZ\\kD_analysis\\R75A_R224A_titrations\\all_CL_reps\\"

files = ["rep1\\AqpZ_R75A_R224A_80CL_unidecfiles", "rep2\\AqpZ_R75A_R224A_80CL_unidecfiles",
         "rep3\\AqpZ_R75A_R224A_80CL_unidecfiles"]

topdir = "Z:\\Group Share\\Megan E\\AqpZ\\kD_analysis\\R75A_R224A_titrations\\all_CL_reps\\"

files = ["rep1\\AqpZ_R75A_R224A_40CL_unidecfiles", "rep2\\AqpZ_R75A_R224A_40CL_unidecfiles",
         "rep3\\AqpZ_R75A_R224A_40CL_unidecfiles"]'''

kds = []
wts = []
mts = []
maxl=5
fig, ax = plt.subplots()
for f in files:
    path = os.path.join(topdir, f, "Extract_grid_2D_Extract.txt")
    print(path, f)
    data = np.loadtxt(path)

    bm = data[:, 1] == 0
    mutant = data[bm][:, 2]
    wt = data[np.logical_not(bm)][:, 2]

    wt = wt[:maxl]
    mutant = mutant[:maxl]
    wt =wt/np.sum(wt)
    mutant = mutant/np.sum(mutant)
    wts.append(wt)
    mts.append(mutant)

    ax.plot(np.arange(len(wt))+0.5, wt, color="g")
    ax.plot(np.arange(len(mutant)), mutant, color="b")

    value = get_kd(mutant, wt, 0, maxl-1)
    kdlist = [value]
    for n in range(0, maxl-1):
        value = get_kd(mutant, wt, n, 1)
        kdlist.append(value)
    kds.append(kdlist)

kds = np.array(kds)
wts = np.array(wts)
mts = np.array(mts)


wavg = np.mean(wts, axis=0)
mavg = np.mean(mts, axis=0)

wstd = np.std(wts, axis=0)
mstd = np.std(mts, axis=0)


print(kds)

maxn=maxl
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

for n in range(1, maxn):
    avgdg = np.mean(kds[:,n,1])
    sddg = np.std(kds[:,n,1])

    text = "\u0394\u0394G = " + str(round(avgdg, 2)) + " \u00B1 " + str(round(sddg, 2))+" "+units
    ax.text(n,wavg[n]+0.2, text)

plt.show()
