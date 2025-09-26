import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from matplotlib import rcParams
from matplotlib import cm
import itertools


path = "F:\\Group Share\\Figures\\Lipid_Tail_comparison\\LL37\\Text Files\\"
os.chdir(path)

#Create empty lists for the names of peptides and the averages that will be plotted
lipidname_Ind3=()
lipidname_Mag3=()
lipidname_TempA3=()
Averages_Ind3=()
Averages_Ind9=()
Averages_Mag3=()
Averages_Mag9=()
Averages_TempA3=()
Averages_TempA9=()
Errors_Ind3=()
Errors_Ind9=()
Errors_Mag3=()
Errors_Mag9=()
Errors_TempA3=()
Errors_TempA9=()

#Load files and append values to empty lists created above
lipidname_Ind3 = np.append(lipidname_Ind3, np.loadtxt("Indolicidin_3to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages_Ind3 = np.append(Averages_Ind3, np.loadtxt("Indolicidin_3to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_Ind3 = np.append(Errors_Ind3, np.loadtxt("Indolicidin_3to1_Errors.txt", delimiter='\t', usecols=(1)))
lipidname_Mag3 = np.append(lipidname_Mag3, np.loadtxt("Magainin_3to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages_Mag3 = np.append(Averages_Mag3, np.loadtxt("Magainin_3to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_Mag3 = np.append(Errors_Mag3, np.loadtxt("Magainin_3to1_Errors.txt", delimiter='\t', usecols=(1)))
lipidname_TempA3 = np.append(lipidname_TempA3, np.loadtxt("TemporinA_3to1_Averages.txt", dtype=str, delimiter='\t', usecols=(0)))
Averages_TempA3 = np.append(Averages_TempA3, np.loadtxt("TemporinA_3to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_TempA3 = np.append(Errors_TempA3, np.loadtxt("TemporinA_3to1_Errors.txt", delimiter='\t', usecols=(1)))
Averages_Ind9 = np.append(Averages_Ind9, np.loadtxt("Indolicidin_9to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_Ind9 = np.append(Errors_Ind9, np.loadtxt("Indolicidin_9to1_Errors.txt", delimiter='\t', usecols=(1)))
Averages_Mag9 = np.append(Averages_Mag9, np.loadtxt("Magainin_9to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_Mag9 = np.append(Errors_Mag9, np.loadtxt("Magainin_9to1_Errors.txt", delimiter='\t', usecols=(1)))
Averages_TempA9 = np.append(Averages_TempA9, np.loadtxt("TemporinA_9to1_Averages.txt", delimiter='\t', usecols=(1)))
Errors_TempA9 = np.append(Errors_TempA9, np.loadtxt("TemporinA_9to1_Errors.txt", delimiter='\t', usecols=(1)))

#Set colormap that will be used in bar chart
cmap = matplotlib.cm.get_cmap('viridis')
colors = cmap(np.arange(0, len(lipidname_Ind3)) / (len(lipidname_Ind3) - 1))

#set the number or subplots in the figure
fig, axs = plt.subplots(3, 2, sharex = 'all')

#Plot the bar charts as subplots
rects1 = axs[0,0].bar(lipidname_Mag3, Averages_Mag3, yerr=Errors_Mag3, color=colors, capsize=3)
rects2 = axs[1,0].bar(lipidname_Ind3, Averages_Ind3, yerr=Errors_Ind3, color=colors, capsize=3)
rects3 = axs[2,0].bar(lipidname_TempA3, Averages_TempA3, yerr=Errors_TempA3, color=colors, capsize=3)
rects4 = axs[0,1].bar(lipidname_Mag3, Averages_Mag9, yerr=Errors_Mag9, color=colors, capsize=3)
rects5 = axs[1,1].bar(lipidname_Ind3, Averages_Ind9, yerr=Errors_Ind9, color=colors, capsize=3)
rects6 = axs[2,1].bar(lipidname_TempA3, Averages_TempA9, yerr=Errors_TempA9, color=colors, capsize=3)

#T-test to determine which lipids have significantly different means using the means and standard deviations
Pairs_Ind9 = list(zip(Averages_Ind9, Errors_Ind9, lipidname_Ind3))
Ind9_combination = list(itertools.combinations(Pairs_Ind9,2))
print('Indolicidin 9:1')
for combo in Ind9_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

Pairs_Mag9 = list(zip(Averages_Mag9, Errors_Mag9,lipidname_Ind3))
Mag9_combination = list(itertools.combinations(Pairs_Mag9, 2))
print('Magainin 9:1')
for combo in Mag9_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

Pairs_TempA9 = list(zip(Averages_TempA9, Errors_TempA9,lipidname_Ind3))
TempA9_combination = list(itertools.combinations(Pairs_Mag9, 2))
print('Temporin A 9:1')
for combo in TempA9_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

Pairs_Ind3 = list(zip(Averages_Ind3, Errors_Ind3, lipidname_Ind3))
Ind3_combination = list(itertools.combinations(Pairs_Ind3, 2))
print('Indolicidin 3:1')
for combo in Ind3_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

Pairs_Mag3 = list(zip(Averages_Mag3, Errors_Mag3, lipidname_Ind3))
Mag3_combination = list(itertools.combinations(Pairs_Mag3, 2))
print('Magainin 3:1')
for combo in Mag3_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

Pairs_TempA3 = list(zip(Averages_TempA3, Errors_TempA3,lipidname_Ind3))
TempA3_combination = list(itertools.combinations(Pairs_Mag3, 2))
print('Temporin A 3:1')
for combo in TempA3_combination:
    t, p = stats.ttest_ind_from_stats(mean1=combo[0][0], std1=combo[0][1], nobs1=3, mean2=combo[1][0], std2=combo[1][1], nobs2=3)
    if p < 0.05 or p == 0.05:
        print(p, combo)

#Place Letter above each bar
Letters = ("A", "B", "C", "D", "E")
for idx,rect in enumerate(rects1):
    height = rect.get_height()
    axs[0,0].text(rect.get_x() + rect.get_width()/2, height + Errors_Mag3[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)

for idx, rect in enumerate(rects2):
    height = rect.get_height()
    axs[1,0].text(rect.get_x() + rect.get_width()/2, height + Errors_Ind3[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)

for idx, rect in enumerate(rects3):
    height = rect.get_height()
    axs[2,0].text(rect.get_x() + rect.get_width()/2, height + Errors_TempA3[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)

for idx,rect in enumerate(rects4):
    height = rect.get_height()
    axs[0,1].text(rect.get_x() + rect.get_width()/2, height + Errors_Mag9[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)

for idx, rect in enumerate(rects5):
    height = rect.get_height()
    axs[1,1].text(rect.get_x() + rect.get_width()/2, height + Errors_Ind9[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)

for idx, rect in enumerate(rects6):
    height = rect.get_height()
    axs[2,1].text(rect.get_x() + rect.get_width()/2, height + Errors_TempA9[idx] + 0.3, Letters[idx], ha='center', rotation=0, size=14)


rcParams['font.size']=40
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42


# Add some text for labels, title and custom x-axis tick labels, etc.
axs[0,0].set_title('3:1', size=26)
#axs[1,0].set_title('Indolicidin (3:1)', size=20)
#axs[2,0].set_title('TemporinA (3:1)', size=20)
axs[0,1].set_title('9:1', size=26)
#axs[1,1].set_title('Indolicidin (9:1)', size=20)
#axs[2,1].set_title('TemporinA (9:1)', size=20)
axs[0,0].set_ylabel('Avg Peptides Incorporated', size=12)
axs[1,0].set_ylabel('Avg Peptides Incorporated', size=12)
axs[2,0].set_ylabel('Avg Peptides Incorporated', size=12)
axs[0,1].yaxis.set_label_position("right")
axs[0,1].yaxis.set_label_coords(-1.42, 0.5)
axs[0,1].set_ylabel('Magainin-2', size=24, rotation=0)
axs[1,1].yaxis.set_label_position("right")
axs[1,1].yaxis.set_label_coords(-1.42, 0.5)
axs[1,1].set_ylabel('Indolicidin', size=24, rotation=0)
axs[2,1].yaxis.set_label_position("right")
axs[2,1].yaxis.set_label_coords(-1.42, 0.5)
axs[2,1].set_ylabel('Temporin A', size=24, rotation=0)

#Set Y-axis range and tick mark spacing
axs[0,0].set_yticks(np.arange(0, 9.05, step=2))
axs[0,0].set_ylim(0, 9.0)
axs[1,0].set_yticks(np.arange(0, 9.05, step=2))
axs[1,0].set_ylim(0, 9.0)
axs[2,0].set_yticks(np.arange(0, 9.05, step=2))
axs[2,0].set_ylim(0, 9.0)
axs[0,1].set_yticks(np.arange(0,9.05, step=2))
axs[0,1].set_ylim(0, 9.0)
axs[1,1].set_yticks(np.arange(0,9.05, step=2))
axs[1,1].set_ylim(0, 9.0)
axs[2,1].set_yticks(np.arange(0,9.05, step=2))
axs[2,1].set_ylim(0,9.0)

#Turn off x-axis labels and tick marks
axs[0,0].get_xaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[2,0].get_xaxis().set_visible(True)
axs[0,1].get_xaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[2,1].get_xaxis().set_visible(True)
axs[0,0].spines['top'].set_visible(False)
axs[1,0].spines['top'].set_visible(False)
axs[2,0].spines['top'].set_visible(False)
axs[0,1].spines['top'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[2,1].spines['top'].set_visible(False)
axs[0,0].spines['right'].set_visible(False)
axs[1,0].spines['right'].set_visible(False)
axs[2,0].spines['right'].set_visible(False)
axs[0,1].spines['right'].set_visible(False)
axs[1,1].spines['right'].set_visible(False)
axs[2,1].spines['right'].set_visible(False)

axs[2,0].tick_params('x', labelsize=20)
axs[2,1].tick_params('x', labelsize=20)
axs[0,0].tick_params('y', labelsize=16)
axs[1,0].tick_params('y', labelsize=16)
axs[2,0].tick_params('y', labelsize=16)
axs[0,1].tick_params('y', labelsize=16)
axs[1,1].tick_params('y', labelsize=16)
axs[2,1].tick_params('y', labelsize=16)
axs[2,0].tick_params(axis='x', length=0)
axs[2,1].tick_params(axis='x', length=0)

plt.show()