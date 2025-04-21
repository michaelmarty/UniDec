import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import scipy.stats as stats

from matplotlib import rcParams
import matplotlib as mpl

mpl.use("TkAgg")

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 11
rcParams['font.sans-serif'] = "Arial"

path = "Z:\\Group Share\\Dash\\Daptomycin\\Oxidized"
os.chdir(path)

file = "FPOPResults.xlsx"
file2 = "NORMALIZED OUTPUT.xlsx"
df = pd.read_excel(file)

equal_var = True

for index, row in df.iterrows():

    pep1 = row["Peptide"]
    m1 = row["Mean"]
    s1 = row["SD"]
    n1 = row["N"]
    for index2, row2 in df.iterrows():
        pep2 = row2["Peptide"]
        newname = row2["Lipid"] + str(row2["Ratio"])
        if pep1 == pep2:
            m2 = row2["Mean"]
            s2 = row2["SD"]
            n2 = row2["N"]
            if n1 > 0 and n2 > 0:
                pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]
            else:
                pval = 1
            if pval < 0.05:
                df.loc[index, newname] = pval
            else:
                df.loc[index, newname] = ""

plt.figure()

colors = ["b", "k", "g", "skyblue", "grey", "limegreen"]
lipids = ["DMPG", "DMPC", "Ecoli"]
peptides = ["Da"]
heights = [0.30, 0.30]

ratios = [3, 9, 18]

ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

for k, p in enumerate(peptides):
    count2 = 0
    for j, l in enumerate(lipids):
        count1 = 0
        for i, r in enumerate(ratios):
            for index, row in df.iterrows():
                if row["Peptide"] == p and row["Lipid"] == l and row["Ratio"] == r:
                    pep1 = row["Peptide"]
                    m1 = row["Mean"]
                    s1 = row["SD"]
                    n1 = row["N"]
                    # cm1 = row["CMean"]
                    # cs1 = row["CSD"]
                    # m1 = m1-cm1
                    # s1 = np.sqrt(s1 ** 2 + cs1 ** 2)
                    lipid = row["Lipid"]
            if p == "LL37":
                axt = ax1
                axb = ax3
            else:
                axt = ax2
                axb = ax4

            color = colors[j]
            color2 = colors[j + 3]

            axt.bar(i + j * 3, m1, yerr=s1, color=color, width=0.4)
            # axt.bar(i + j * 3, cm1, yerr=cs1, color=color2, width=0.4)
            axt.spines['top'].set_visible(False)
            axt.spines['right'].set_visible(False)

            axb.bar(j + i * 3, m1, yerr=s1, color=color, width=0.4)
            # axb.bar(j + i * 3, cm1, yerr=cs1, color=color2, width=0.4)
            axb.spines['top'].set_visible(False)
            axb.spines['right'].set_visible(False)

            h = heights[k]
            for index2, row2 in df.iterrows():
                pep2 = row2["Peptide"]
                tl = row2["Lipid"]
                r2 = row2["Ratio"]
                m2 = row2["Mean"]
                s2 = row2["SD"]
                n2 = row2["N"]

                # cm2 = row["CMean"]
                # cs2 = row["CSD"]
                # m2 = m2 - cm2
                # s2 = np.sqrt(s2**2 + cs2 **2)
                if n1 > 0 and n2 > 0:
                    pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]
                else:
                    pval = 1
                spacing = 0.02
                if pep1 == pep2 and tl == lipid and r2 > r and pval < 0.05:
                    if r2 - r == 6 or r2 - r == 9:
                        w = 1
                    else:
                        w = 2
                    x1 = i + j * 3 + spacing * count1
                    x2 = x1 + w
                    h = h - 0.005 * count1
                    axt.vlines(x1, m1 + s1 + 0.01, h)
                    axt.vlines(x2, m2 + s2 + 0.01, h)
                    axt.hlines(h, x1, x2)
                    count1 = count1 + 1

                for li, lip in enumerate(lipids):
                    if tl == lip:
                        lindex = li
                if pep1 == pep2 and r2 == r and lindex > j and pval < 0.05:
                    w = lindex - j
                    x1 = j + i * 3 + spacing * count2
                    x2 = x1 + w
                    h = h - 0.003 * count2
                    axb.vlines(x1, m1 + s1 + 0.01, h)
                    axb.vlines(x2, m2 + s2 + 0.01, h)
                    axb.hlines(h, x1, x2)
                    count2 = count2 + 1
                    print(count2, lindex, j, tl)

plt.savefig("OUTPUT.pdf")
#plt.savefig("OUTPUT.png")
plt.show()

#df.to_excel(file2, index=False)
