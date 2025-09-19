
import matplotlib.pyplot as plt
import pandas as pd
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

file =
file2 =
df = pd.read_excel(file)


equal_var = True

for index, row in df.iterrows():
    pep1 = row["Natural Lipids"]
    m1 = row["Mean"]
    s1 = row["SD"]
    n1 = row["N"]

plt.figure()

colors = ["y", "b", "k", "g", "y", "skyblue", "grey", "limegreen"]
sizes = ["MSP1D1", "MSP1E3D1" ]
NaturalLipids = ["Ecoli"]
heights = [45]
hspacing = 0.05
temperatures = [4, 25, 41]

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)


for k, p in enumerate(NaturalLipids):
    count2 = 0
    for j, l in enumerate(sizes):
        count1 = 0
        for i, r in enumerate(temps):
            for index, row in df.iterrows():
                if row["Natural Lipids"] == p and row["Lipid"] == l and row["Ratio"] == r:
                    pep1 = row["Natural Lipids"]
                    m1 = row["Mean"]
                    s1 = row["SD"]
                    n1 = row["N"]

                    size = row["Size"]

            axt = ax1
            axb = ax2

            color = colors[j]
            #color2 = colors[j + 4]

            axt.bar(i + 0.5 + j * len(temps), m1, yerr=s1, color=color, width=0.4)
            #axt.bar(i + j * 3, cm1, yerr=cs1, color=color2, width=0.4)
            axt.spines['top'].set_visible(False)
            axt.spines['right'].set_visible(False)

            axb.bar(j + i * len(sizes) + 0.5, m1, yerr=s1, color=color, width=0.4)
            #axb.bar(j + i * 3, cm1, yerr=cs1, color=color2, width=0.4)
            axb.spines['top'].set_visible(False)
            axb.spines['right'].set_visible(False)

            toph = heights[k]
            for index2, row2 in df.iterrows():
                pep2 = row2["Natural Lipids"]
                tl = row2["Size"]
                r2 = row2["Temperature"]
                m2 = row2["Mean"]
                s2 = row2["SD"]
                n2 = row2["N"]
                pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]

                if pep1 == pep2 and tl == lipid and r2 > r and pval<0.05:
                    if r2 - r == 6 or r2 - r == 9:
                        w = 1
                    else:
                        w = 2
                    x1 = i + 0.5 + j * 3 + hspacing * count1
                    x2 = x1 + w
                    h = toph - 0.03 * toph * count1
                    axt.vlines(x1, m1 + s1 + 0.01*toph, h)
                    axt.vlines(x2, m2 + s2 + 0.01*toph, h)
                    axt.hlines(h, x1, x2)
                    count1 = count1 + 1

            toph = heights[k]
            for index2, row2 in df.iterrows():
                pep2 = row2["Natural Lipids"]
                tl = row2["Size"]
                r2 = row2["Temperature"]
                m2 = row2["Mean"]
                s2 = row2["SD"]
                n2 = row2["N"]
                pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]
                print(tl, l, r, r2, pval, m1, s1, n1, m2, s2, n2)
                for li, lip in enumerate(sizes):
                    if tl == lip:
                        lindex = li
                if pep1 == pep2 and r2 == r and lindex > j and pval < 0.05:
                    w = lindex - j
                    x1 = j + 0.5 + i * len(sizes) + hspacing * (count2%3-1)
                    x2 = x1 + w
                    h = toph - 0.01 * toph * (count2%5-2)
                    axb.vlines(x1, m1 + s1 + 0.01*toph, h)
                    axb.vlines(x2, m2 + s2 + 0.01*toph, h)
                    axb.hlines(h, x1, x2)
                    count2 = count2 + 1
                    print(count2, lindex, j, tl)

#plt.savefig("Z:\\Group Share\\Deseree\\POP Paper 2\\FPOP2 SuperT Output.pdf")
plt.show()

#df.to_excel(file2, index=False)
