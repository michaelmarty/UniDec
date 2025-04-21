import os
import numpy as np
import matplotlib.pyplot as plt

# Define the top-level directory
topdir = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\figure"
file1 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\20240409 TD Native SEC Peak 1 _dilute_unidecfiles\\20240409 TD Native SEC Peak 1 _dilute_input.dat"
file2 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\figure\\output.txt"

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)

x1, y1 = data1[:, 0], (data1[:, 1]/np.amax(data1[:, 1])*100)
x2, y2 = data2[:, 0], data2[:, 1]
plt.plot(x1, y1, label="Ensemble nMS", color='green', alpha=0.5, linewidth=0.5)
plt.plot(x2, y2, label="CDMS", color='blue', alpha=1)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("m/z", fontstyle='italic')
plt.ylabel("%")
plt.xlim(2473, 8000)
plt.ylim(0, 100)
plt.tight_layout()
plt.legend()
plt.savefig(os.path.join(topdir, "output.png"), dpi=300)
plt.show()
