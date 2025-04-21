import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats

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



plt.figure()
x = ["200", "100","50"]
y = [0.40, 1.40, 1.11]

plt.bar(x, y)
plt.gca().set_ylim(-0.5, 2.5)
plt.yticks([-0.5, 0, 0.5, 1.0, 1.5, 2, 2.5])
e = [0.09, 0.07, 0.09]
plt.errorbar(x, y, yerr=e, color='Black')
plt.xlabel("Number of Lipids")
plt.ylabel("Delta G")
plt.savefig(os.path.join("Z:\\Group Share\\Hiruni Jayasekera\\HSJ_FAM_MAE_MutantCyclePaper\\Revisions", "VH_Plot2.pdf"))
plt.show()

