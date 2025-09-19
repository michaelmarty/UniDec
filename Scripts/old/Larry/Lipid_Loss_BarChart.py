import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "C:\\Users\\ldogg\\OneDrive\\Documents\\Marty Lab\\Papers\\"
os.chdir(path)

#file1 = "Delta_Lipid.txt"
#file2 = "Delta_Lipid_Errors.txt"
ydata1 = np.loadtxt("Delta_Lipid.txt")
yerr1 = np.loadtxt("Delta_Lipid_Errors.txt")

fig, ax = plt.subplots()

ind = np.arange(9)
width = 0.5
p1 = ax.bar(ind, ydata1, width, bottom=0, yerr=yerr1, capsize=3.0)

labels = ("Bactenecin", "Reduced Bactenecin", "Indolicidin", "KR12", "Lactoferricin B", "Reduced Lactoferricin B", "Magainin2", "Temporin A", "Temporin L")
x = np.arange(len(labels))
#rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Δ Lipids', size=24)
#ax.set_xlabel('Peptide', size=24)
ax.set_yticklabels(np.arange(0, 41, step=(10)))
#ax.set_title('DMPC', size=24)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45)
#ax.legend()
plt.xticks(size=10)
plt.xlim(-0.5,len(labels))
plt.ylim(top=45)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,41, step=10), size=20)
ax.autoscale_view()
fig.tight_layout()

plt.show()