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

file = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\Flash Extraction MP Nanodisc\\Data\\AqpZ\\Fig5Data.xlsx"
file = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\Flash Extraction MP Nanodisc\\Data\\AqpZ\\JK data AmtB Fig 1 and Fig 3.xlsx"

df = pd.read_excel(file)

lipids = np.unique(df["Lipid Combo"])

print(df)
print(lipids)

equal_var = True

for index, row in df.iterrows():
    l1 = row["Lipid Combo"]
    c1 = row["Condition"]
    m1 = row["Mall"]
    s1 = row["Sall2"]
    n1 = row["N"]*3

    for index2, row2 in df.iterrows():
        c2 = row2["Condition"]
        l2 = row2["Lipid Combo"]
        if c1 != c2 and l1==l2:
            m2 = row2["Mall"]
            s2 = row2["Sall2"]
            n2 = row2["N"]*3
            pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]
            if pval < 0.05:
                print(l1, c1, c2, pval)
                pass
            else:
                #print(l1, c1, c2, pval, "Not Significant")
                pass




#pval = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=equal_var)[1]