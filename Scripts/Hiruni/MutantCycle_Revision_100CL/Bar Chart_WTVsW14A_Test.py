import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats
import pandas as pd

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

df = pd.read_excel("Z:\\Group Share\\Hiruni Jayasekera\\HSJ_FAM_MAE_MutantCyclePaper\\MCPaper_Revisions.xlsx", sheet_name='Lipid_1')

print(df)

labels = ["CE_L1", "W14A_R224A_L1"]

figure, subplots = plt.subplots(1, len(labels))

for i, l in enumerate(labels):
    ax = subplots[i]
    y = df[l]
    if i == 0:
        e = df["E_L1"]
    else:
        e = df["E_L1." + str(i)]

    ax.set_ylim(-0.1, 1.8)
    ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0, 1.5])
    if i!=0:
        ax.set_yticklabels(["","","","","",""])
    ax.bar([0, 1, 2], y, color = ["#C1FFC1", "#32CD32", "#006400"])
    ax.errorbar([0, 1, 2], y, yerr=e, color='Black', linestyle=" ")
    ax.set_xticks([0, 1, 2], ["50", "100", "200"])

    ax.set_xlabel("CL Ratio")
    if i ==0:
        ax.set_ylabel("ΔΔG (kJ/mol)")
    ax.set_title(l)

plt.tight_layout()
plt.savefig(os.path.join("Z:\\Group Share\\Hiruni Jayasekera\\HSJ_FAM_MAE_MutantCyclePaper\\Revisions", "L1_CE_W14A_R224A.pdf"))
plt.show()

