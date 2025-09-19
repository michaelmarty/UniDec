import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
from matplotlib import cm


path = "D:\\Group Share\\Figures\\Peptide_Paper2\\All Peptides\\Weighted_Average\\"
os.chdir(path)

#Create empty lists for the names of peptides and the averages that will be plotted
Pepname3=()
Pepname9=()
Pepname18=()
Averages3=()
Averages9=()
Averages18=()
Errors3=()
Errors9=()
Errors18=()

#Load files and append values to empty list created above
Pepname3 = np.append(Pepname3, np.loadtxt("3to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages3 = np.append(Averages3, np.loadtxt("3to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors3 = np.append(Errors3, np.loadtxt("3to1_Errors.txt", delimiter='\t', usecols=(1)))
Pepname9 = np.append(Pepname9, np.loadtxt("9to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages9 = np.append(Averages9, np.loadtxt("9to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors9 = np.append(Errors9, np.loadtxt("9to1_Errors.txt", delimiter='\t', usecols=(1)))
Pepname18 = np.append(Pepname18, np.loadtxt("18to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages18 = np.append(Averages18, np.loadtxt("18to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors18 = np.append(Errors18, np.loadtxt("18to1_Errors.txt", delimiter='\t', usecols=(1)))

#Set colormap that will be used in bar chart
cmap = matplotlib.cm.get_cmap('viridis')
colors = cmap(np.arange(0, len(Pepname3)) / (len(Pepname3) - 1))
colors9 = cmap(np.arange(0, len(Pepname9)) / (len(Pepname9) - 1))
colors18 = cmap(np.arange(0, len(Pepname18)) / (len(Pepname18) - 1))

#set the number or subplots in the figure
fig, (ax1, ax2, ax3) = plt.subplots(3)
#fig, (ax1, ax2) = plt.subplots(2)

#Plot the bar charts as subplots
rects1 = ax1.bar(Pepname3, Averages3, yerr=Errors3, color=colors)
rects2 = ax2.bar(Pepname9, Averages9, yerr=Errors9, color=colors9)
rects3 = ax3.bar(Pepname18, Averages18, yerr=Errors18, color=colors18)

#Place Peptide name above each bar
for idx,rect in enumerate(rects1):
    height = rect.get_height()
    ax1.text(rect.get_x() + rect.get_width()/2., height, Pepname3[idx], ha='center', va='bottom', rotation=90, size=14)

for idx, rect in enumerate(rects2):
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width() / 2., height, Pepname9[idx], ha='center', va='bottom', rotation=90, size=14)

for idx, rect in enumerate(rects3):
    height = rect.get_height()
    ax3.text(rect.get_x() + rect.get_width() / 2., height, Pepname18[idx], ha='center', va='bottom', rotation=90, size=14)

#rcParams['font.size']=20
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

# Add some text for labels, title and custom x-axis tick labels, etc.
ax1.set_ylabel('Avg Peptides Incorporated', size=12)
ax2.set_ylabel('Avg Peptides Incorporated', size=12)
ax3.set_ylabel('Avg Peptides Incorporated', size=12)
ax1.set_title('Low (3:1)', size=20)
ax2.set_title('Intermediate (9:1)', size=20)
ax3.set_title('High (18:1*)', size=20)
#ax1.set_xticklabels(Pepname3, rotation=90, size=20)
#ax2.set_xticklabels(Pepname9, rotation=90, size=20)
#ax3.set_xticklabels(Pepname9, rotation=90, size=20)

#Set Y-axis range and tick mark spacing
ax1.set_yticks(np.arange(0, 8.6, step=2))
ax1.set_ylim(0, 8.5)
ax2.set_yticks(np.arange(0, 8.6, step=2))
ax2.set_ylim(0, 8.5)
ax3.set_yticks(np.arange(0, 8.6, step=2))
ax3.set_ylim(0, 8.5)

#Turn off x-axis labels and tick marks
ax1.tick_params(axis='x', labelbottom=False)
ax2.tick_params(axis='x', labelbottom=False)
ax3.tick_params(axis='x', labelbottom=False)
ax1.get_xaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax3.get_xaxis().set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax3.spines['right'].set_visible(False)

#fig.tight_layout()

plt.show()