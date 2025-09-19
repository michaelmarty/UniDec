import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
from matplotlib import cm
import seaborn as sb


path = "F:\\Group Share\\Figures\\Peptide_Paper2\\All Peptides\\Weighted_Average\\"
os.chdir(path)
Pepname3=()
Pepname9=()
Averages3=()
Averages9=()
Errors3=()
Errors9=()
Averages18=()
Errors18=()
Pepname3 = np.append(Pepname3, np.loadtxt("3to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages3 = np.append(Averages3, np.loadtxt("3to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors3 = np.append(Errors3, np.loadtxt("3to1_Errors.txt", delimiter='\t', usecols=(1)))
Pepname9 = np.append(Pepname9, np.loadtxt("9to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages9 = np.append(Averages9, np.loadtxt("9to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors9 = np.append(Errors9, np.loadtxt("9to1_Errors.txt", delimiter='\t', usecols=(1)))
Averages18 = np.append(Averages18, np.loadtxt("18to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors18 = np.append(Errors18, np.loadtxt("18to1_Errors.txt", delimiter='\t', usecols=(1)))

#cmap = matplotlib.cm.get_cmap('viridis')
#colors = cmap(np.arange(0, len(Pepname3)) / (len(Pepname3) - 1))
#colors9 = cmap(np.arange(0, len(Pepname9)) / (len(Pepname9) - 1))

fig, (ax1, ax2, ax3) = plt.subplots(3)

ax1.sb.barplot(Pepname3, Averages3)
ax2.sb.barplot(Pepname9, Averages9)
ax3.sb.barplot(Pepname9, Averages18)
#ax1.bar(Pepname3, Averages3, yerr=Errors3, color=colors)
#ax2.bar(Pepname9, Averages9, yerr=Errors9, color=colors9)
#ax3.bar(Pepname9, Averages18, yerr=Errors18, color=colors9)




rcParams['font.size']=20
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

# Add some text for labels, title and custom x-axis tick labels, etc.
#ax1.set_ylabel('AVG Peptides Incorporated')
#ax2.set_ylabel('AVG Peptides Incorporated')
#ax3.set_ylabel('AVG Peptides Incorporated')
#ax1.set_title('Low (3:1)')
#ax2.set_title('Intermediate (9:1)')
#ax3.set_title('High (18:1*)')
#ax1.set_xticks(Pepname3)
#ax2,ax3.set_xticks(Pepname9)
#ax1.set_yticks(np.arange(0, 2.6, step=0.5))
#ax1.set_ylim(0, 3)
#ax2.set_yticks(np.arange(0, 7.1, step=1))
#ax2.set_ylim(0, 7)
#ax3.set_yticks(np.arange(0, 9.1, step=1.5))
#ax3.set_ylim(0, 9)
#plt.yticks(size=20)
#plt.xticks(size=20)
#plt.ylim(top=6.5)
#plt.ylim(bottom=0.0)
#fig.tight_layout()
#fig2.tight_layout()

plt.show()