import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl

mpl.use("WxAgg")

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

#input file containing folder. input file must be in .csv format.
path = "Z:\\Group Share\\Dash\\Daptomycin\\Oxidized"
os.chdir(path)
file = "Output.csv"

df = pd.read_csv(file)
df.columns = ["Path", "Lipid", "Ratio", "Ox"]

lipids = ["DMPC", "DMPG", "Ecoli", "Free", "DPPC", "DLPC", "POPC", "fdap", "fdapb"]
ratios = [3, 9, 18]

colors = ["orange", "b", "g", "y", "r", "purple", "pink", "k", "grey", "k"]

output = []
for i, l in enumerate(lipids):
    for j, r in enumerate(ratios):
        df2 = df[df["Lipid"]==l]
        df2 = df2[df2["Ratio"]==r]

        ox = df2["Ox"].to_numpy()
        n = len(ox)
        if n > 0:
            m = np.mean(ox)
            #s = np.std(ox, ddof=1)
            out = ["Da", l, r, m, s, n]
        else:
            out = ["Da", l, r, 0, 0, n]
        output.append(out)
        plt.bar(i*3+j, m, yerr=s, color=colors[i])


df3 = pd.DataFrame(output)
df3.columns = ["Peptide", "Lipid", "Ratio", "Mean", "SD", "N"]
print(df3)
df3.to_excel("FPOPResults.xlsx")
plt.ylim(0,0.8)
plt.show()