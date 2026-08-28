import matplotlib as mpl
mpl.use("Agg")

import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.size'] = 12

dir = r"Z:\Group Share\BHT\Temp\Paper\DMT\Protein Mix 21 Injections Set 1\merged\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_unidecfiles"
os.chdir(dir)

kfile = "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_htkernel.txt"
kdat = np.loadtxt(kfile)

tfile = "TIC_chrom.txt"
tdat = np.loadtxt(tfile)

hfile = "TIC_decon.txt"
hdat = np.loadtxt(hfile)

padindex = 0
rollindex = 0
ymax = 900  # y-axis top value

for i in range(0, 588):
    startindex = i

    # Shift kernel and pad
    kdat2 = np.concatenate([np.zeros(padindex), np.roll(kdat, startindex)])

    # Scale red line to reach top of ylim
    if np.max(kdat2) > 0:
        kdat2 = kdat2 / np.max(kdat2) * ymax

    print(i)
    plt.figure(figsize=(16,9))
    plt.plot(kdat2, color="red")
    plt.plot(tdat[:,1], color="black")

    hindex = startindex + padindex + rollindex
    plt.plot(hdat[:hindex,1], color="blue")

    plt.xlabel("Scan Number")
    plt.ylabel("Ion Count")
    plt.ylim(0, ymax)

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig("C:\\Python\\Animations\\An" + str(startindex) + ".png",
                dpi=300, bbox_inches="tight")
    plt.close()
