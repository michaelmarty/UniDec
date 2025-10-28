import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
#rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=1
rcParams['axes.linewidth']=1
rcParams['font.size']=12

x = [0,
     1.5,
     3,
     6,
     9,
     12,
     18,
     24,
     42,
     60,
     80
     ]
'''
y = [
    0.951,
    0.918,
    0.885,
    0.681,
    0.450,
    0.321,
    0.204,
    0.220,
    0.139,
    0.222,
    0.190
]
'''
yerr = [
    0.034695837,
    0,
    0.096697941,
    0.073189531,
    0.076315568,
    0.045549811,
    0.062030047,
    0.085524175,
    0,
    0,
    0
]

files=["Z:\\Group Share\\Larry Walker\\MassSpecData\\Melittin\\20181011_metittinDMPG\\UniDec_Figures_and_Files\\MelittinTitraion_DMPG_1010_Mass_Defect_Extracts.txt"]
for f in files:
    data = np.loadtxt(f)
    print(data)

for i in range(0,len(data)):
    x=np.arange(0,len(data[i]))
    y=np.array(data[i])
    yerr=y*0.1

    plt.errorbar(x, y, yerr=yerr, marker="o")
#plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
#plt.errorbar(x, y-0.4, yerr=yerr, marker="o")

plt.show()
