import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.5
rcParams['axes.linewidth'] = 0.5
# rcParams['font.size']=8


R = 8.31446261815324 / 1000.
units = "kJ / mol"
T = 303.15
R = 1.98720425864083 / 1000.
units = "kcal / mol"

file = "Z:\Group Share\Megan E\AqpZ\kD_analysis\R75A_titrations\\20211012_POPG_RT\\20211012_MAE_AqpZ_R75A_WT_80_POPG_01_unidecfiles\\Extract_grid_2D_Extract.txt"
# file = "Z:\\Group Share\\Megan E\\LipidBindingData\\20210830_MAE_AqpZ_WT_R75A_CL_442_unidecfiles\\Extract_grid_2D_Extract.txt"
file = "Z:\Group Share\Megan E\AqpZ\kD_analysis\R224A_titrations\\all_CL_reps\\rep3\\20211020_mae_AqpZ_R224A_80CL_2_02_unidecfiles\\Extract_grid_2D_Extract.txt"

data = np.loadtxt(file)

bm = data[:, 1] == 0
mutant = data[bm][:, 2]

wt = data[np.logical_not(bm)][:, 2]


def get_kd(mutant, wt, n=0, l=1):
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    deltag = -R * T * np.log(kd)
    return kd, deltag, deltag / l


fig, ax = plt.subplots()
print("Total", get_kd(mutant, wt, 0, 4), units)
for n in range(0, 4):
    value = get_kd(mutant, wt, n, 1)
    ax.bar(n, mutant[n] / np.amax(mutant), 0.45, color="b")
    ax.bar(n + 0.45, wt[n] / np.amax(wt), 0.45, color="g")
    print(n, " =", value, units)
ax.set_xticks(np.arange(0, 4) + 0.225, labels=["1", "2", "3", "4"])
ax.set_yticks([0, 0.5, 1], labels=["0", "%", "100"])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylim(0,1)
ax.set_xlabel("Number of Bound Lipids")
plt.show()
