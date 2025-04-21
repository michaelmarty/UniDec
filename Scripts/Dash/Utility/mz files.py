import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 18
rcParams['font.weight'] = 500
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 500
rcParams['axes.titlesize'] = 22
rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

file1 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_unidecfiles\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_input.dat"
file1 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_unidecfiles\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_input.dat"
x = np.unique(np.loadtxt(file1,usecols=(0,)))
y = np.unique(np.loadtxt(file1,usecols=(1,)))
z = np.reshape(np.loadtxt(file1,usecols=(2,)), (len(x), len(y)))

z = np.sum(z, axis=1)

znorm =  z / np.max(z)*100

output_file = os.path.join(os.path.dirname(file1), "output.txt")  # Specify the output file path
np.savetxt(output_file, np.column_stack((x, znorm)), fmt='%g', comments='')

plt.plot(x, znorm)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("m/z", fontstyle='italic')
plt.ylabel("Intensity (%)")
plt.xlim(np.min(x), 8000)
plt.ylim(0, 100)
plt.tight_layout()
plt.savefig(os.path.join((os.path.dirname(file1)), "output.png"), dpi=300)
plt.show()




