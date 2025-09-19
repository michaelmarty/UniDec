import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import scipy.stats as stats

from matplotlib import rcParams
import matplotlib as mpl
import TukeyPytest as TPT

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

file = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\ForSuperTTesterRaw_June28.xlsx"
file2 = "Z:\\Group Share\\Dash\\Raw_SuperTukeyTesterOutput.xlsx"
df = pd.read_excel(file)

if True:
    tdf2 = pd.read_csv("Z:\\Group Share\\Deseree\\FPOP Paper 2\\4R_Raw_June14_Tukey.csv")
    print(tdf2)
if False:
    tdf = pd.read_csv("Z:\\Group Share\\Deseree\\FPOP Paper 2\\4R_Raw_June14.csv")
    tdf2 = TPT.run_tukey(tdf)
    print(tdf2)

equal_var = True

for index, row in df.iterrows():
    pep1 = row["Peptide"]
    m1 = row["Mean"]
    s1 = row["SD"]
    n1 = row["N"]

plt.figure()

colors = ["cyan", "b", "k", "g", "purple", "orange", "c", "yellow"]
lipids = ["DLPC", "DMPG", "DMPC", "E.coli", "DPPC", "POPC", "DLPC", "None"]
peptides = ["Daptomycin"]
heights = [70]
hspacing = 0.05
ratios = [3, 9, 18]

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

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

            axt = ax1
            axb = ax2

            color = colors[j]
            axt.bar(i + 0.5 + j * len(ratios), m1, yerr=s1, color=color, width=0.4)
            axt.spines['top'].set_visible(False)
            axt.spines['right'].set_visible(False)

            axb.bar(j + i * len(lipids) + 0.5, m1, yerr=s1, color=color, width=0.4)
            axb.spines['top'].set_visible(False)
            axb.spines['right'].set_visible(False)

            toph = heights[k]
            for index2, row2 in df.iterrows():
                pep2 = row2["Peptide"]
                tl = row2["Lipid"]
                r2 = row2["Ratio"]
                m2 = row2["Mean"]
                s2 = row2["SD"]

                testval = TPT.check_tukey(tdf2, l, r, tl, r2)

                if pep1 == pep2 and tl == l and r2 > r and testval:
                    if r2 - r == 6 or r2 - r == 9:
                        w = 1
                    else:
                        w = 2
                    x1 = i + 0.5 + j * 3 + hspacing * count1
                    x2 = x1 + w
                    h = toph - 0.03 * toph * count1
                    axt.vlines(x1, m1 + s1 + 0.01 * toph, h)
                    axt.vlines(x2, m2 + s2 + 0.01 * toph, h)
                    axt.hlines(h, x1, x2)
                    count1 = count1 + 1

                for li, lip in enumerate(lipids):
                    if tl == lip:
                        lindex = li
                if pep1 == pep2 and r2 == r and lindex > j and testval:
                    w = lindex - j
                    x1 = j + 0.5 + i * len(lipids) + hspacing * (count2 % 3 - 1)
                    x2 = x1 + w
                    h = toph - 0.01 * toph * (count2 % 5 - 2)
                    axb.vlines(x1, m1 + s1 + 0.01 * toph, h)
                    axb.vlines(x2, m2 + s2 + 0.01 * toph, h)
                    axb.hlines(h, x1, x2)
                    count2 = count2 + 1
                    print(count2, lindex, j, tl)

plt.savefig("Z:\\Group Share\\Deseree\\FPOP Paper 2\\NEW_SuperTukeyTesterOutput.pdf")
plt.savefig("Z:\\Group Share\\Deseree\\FPOP Paper 2\\NEW_SuperTukeyTesterOutput.png")
plt.show()

# df.to_excel(file2, index=False)
