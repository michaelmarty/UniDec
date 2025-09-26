import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from matplotlib import rcParams
from matplotlib import cm
import itertools


path = "Z:\Julia\M2 protein\wildtype/NEW PURIFICATION/drugbinding/c8e4 9/text files"
os.chdir(path)



#Create empty lists for the names of peptides and the averages that will be plotted
nodrug=()
nineteen=()
threeseven=()
sevenfive=()
hundredfifty=()

c8e4_9_0drug=()
c8e4_9_18drug=()
c8e4_9_37drug=()
c8e4_9_75drug=()
c8e4_9_150drug=()
s31n_c8e4_9_0drug=()
s31n_c8e4_9_18drug=()
s31n_c8e4_9_37drug=()
s31n_c8e4_9_75drug=()
s31n_c8e4_9_150drug=()

error_c8e4_9_0drug=()
error_c8e4_9_18drug=()
error_c8e4_9_37drug=()
error_c8e4_9_75drug=()
error_c8e4_9_150drug=()
error_s31n_c8e4_9_0drug=()
error_s31n_c8e4_9_18drug=()
error_s31n_c8e4_9_37drug=()
error_s31n_c8e4_9_75drug=()
error_s31n_c8e4_9_150drug=()
#Load files and append values to empty lists created above


nodrug = np.append(c8e4_9_0drug, np.loadtxt("c8e4_9_0drug.txt", dtype=str, delimiter='\t'))
c8e4_9_0drug = np.append(c8e4_9_0drug, np.loadtxt("c8e4_9_0drug.txt", delimiter='\t'))
error_c8e4_9_0drug = np.append(error_s31n_c8e4_9_0drug, np.loadtxt("error_c8e4_9_0drug.txt", delimiter='\t'))
nineteen = np.append(c8e4_9_18drug, np.loadtxt("c8e4_9_18drug.txt", dtype=str, delimiter='\t'))
c8e4_9_18drug = np.append(c8e4_9_18drug, np.loadtxt("c8e4_9_18drug.txt", delimiter='\t'))
error_c8e4_9_18drug = np.append(error_c8e4_9_18drug, np.loadtxt("error_c8e4_9_18drug.txt", delimiter='\t'))
threeseven = np.append(c8e4_9_37drug, np.loadtxt("c8e4_9_37drug.txt", dtype=str, delimiter='\t'))
c8e4_9_37drug = np.append(c8e4_9_37drug, np.loadtxt("c8e4_9_37drug.txt", delimiter='\t'))
error_c8e4_9_37drug = np.append(error_c8e4_9_37drug, np.loadtxt("error_c8e4_9_37drug.txt", delimiter='\t'))
sevenfive = np.append(c8e4_9_75drug, np.loadtxt("c8e4_9_75drug.txt", dtype=str, delimiter='\t'))
c8e4_9_75drug = np.append(c8e4_9_75drug, np.loadtxt("c8e4_9_75drug.txt", delimiter='\t'))
error_c8e4_9_75drug = np.append(error_c8e4_9_75drug, np.loadtxt("error_c8e4_9_75drug.txt", delimiter='\t'))
hundredfifty = np.append(c8e4_9_150drug, np.loadtxt("c8e4_9_150drug.txt", dtype=str, delimiter='\t'))
c8e4_9_150drug = np.append(c8e4_9_150drug, np.loadtxt("c8e4_9_150drug.txt", delimiter='\t'))
error_c8e4_9_150drug = np.append(error_c8e4_9_150drug, np.loadtxt("error_c8e4_9_150drug.txt", delimiter='\t'))



s31n_c8e4_9_0drug = np.append(s31n_c8e4_9_0drug, np.loadtxt("s31n_c8e4_9_0drug.txt", delimiter='\t'))
error_s31n_c8e4_9_0drug = np.append(error_s31n_c8e4_9_0drug, np.loadtxt("error_s31n_c8e4_9_0drug.txt", delimiter='\t'))

s31n_c8e4_9_18drug = np.append(s31n_c8e4_9_18drug, np.loadtxt("s31n_c8e4_9_18drug.txt", delimiter='\t'))
error_s31n_c8e4_9_18drug = np.append(error_s31n_c8e4_9_18drug, np.loadtxt("error_s31n_c8e4_9_18drug.txt", delimiter='\t'))

s31n_c8e4_9_37drug = np.append(s31n_c8e4_9_37drug, np.loadtxt("s31n_c8e4_9_37drug.txt", delimiter='\t'))
error_s31n_c8e4_9_37drug = np.append(error_s31n_c8e4_9_37drug, np.loadtxt("error_s31n_c8e4_9_37drug.txt", delimiter='\t'))

s31n_c8e4_9_75drug = np.append(s31n_c8e4_9_75drug, np.loadtxt("s31n_c8e4_9_75drug.txt", delimiter='\t'))
error_s31n_c8e4_9_75drug = np.append(error_s31n_c8e4_9_75drug, np.loadtxt("error_s31n_c8e4_9_75drug.txt", delimiter='\t'))

s31n_c8e4_9_150drug = np.append(s31n_c8e4_9_150drug, np.loadtxt("s31n_c8e4_9_150drug.txt", delimiter='\t'))
error_s31n_c8e4_9_150drug = np.append(error_s31n_c8e4_9_150drug, np.loadtxt("error_s31n_c8e4_9_150drug.txt", delimiter='\t'))






#Set colormap that will be used in bar chart

cmap = matplotlib.cm.get_cmap('viridis_r')
colors = cmap(np.arange(0, len(nodrug)) / (len(nodrug) - 1))

#plt.figure()
#barlist=plt.bar([1,2,3,4,5,6], [1,2,3,4,5,6])
#barlist[0].setcolor('r')
#barlist[1].setcolor('y')
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1),
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1))

#plt.bar = colors = ['black', 'red', 'green', 'blue', 'cyan', 'magenta']




#set the number or subplots in the figure
fig, axs = plt.subplots(2, 5, sharex = 'all')

#Plot the bar charts as subplots
rects1 = axs[0,0].bar(nodrug, c8e4_9_0drug, yerr= error_c8e4_9_0drug, capsize=3, color=colors)
rects2 = axs[0,1].bar(nineteen, c8e4_9_18drug, yerr= error_c8e4_9_18drug, capsize=3, color=colors)
rects3 = axs[0,2].bar(threeseven, c8e4_9_37drug, yerr= error_c8e4_9_37drug, capsize=3, color=colors)
rects4 = axs[0,3].bar(sevenfive, c8e4_9_75drug, yerr= error_c8e4_9_75drug, capsize=3, color=colors)
rects5 = axs[0,4].bar(hundredfifty, c8e4_9_150drug, yerr= error_c8e4_9_150drug, capsize=3, color=colors)
rects6 = axs[1,0].bar(nodrug, s31n_c8e4_9_0drug, yerr= error_s31n_c8e4_9_0drug, capsize=3, color=colors)
rects7 = axs[1,1].bar(nineteen, s31n_c8e4_9_18drug, yerr= error_c8e4_9_18drug, capsize=3, color=colors)
rects8 = axs[1,2].bar(threeseven, s31n_c8e4_9_37drug, yerr= error_s31n_c8e4_9_37drug, capsize=3, color=colors)
rects9 = axs[1,3].bar(sevenfive, s31n_c8e4_9_75drug, yerr= error_s31n_c8e4_9_75drug, capsize=3, color=colors)
rects10 = axs[1,4].bar(hundredfifty, s31n_c8e4_9_150drug, yerr= error_s31n_c8e4_9_150drug, capsize=3, color=colors)






axs[0, 0].get_yaxis().set_ticks([0, 50, 100])
axs[0, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20)
axs[0,1].get_yaxis().set_ticks([0, 50, 100])
axs[0,2].get_yaxis().set_ticks([0, 50, 100])
axs[0,3].get_yaxis().set_ticks([0, 50, 100])
axs[0,4].get_yaxis().set_ticks([0, 50, 100])

axs[1, 0].get_yaxis().set_ticks([0, 50, 100])
axs[1, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[1,1].get_yaxis().set_ticks([0, 50, 100])
axs[1,2].get_yaxis().set_ticks([0, 50, 100])
axs[1,3].get_yaxis().set_ticks([0, 50, 100])
axs[1,4].get_yaxis().set_ticks([0, 50, 100])



axs[0,1].get_yaxis().set_ticklabels(["", '', ""]),
axs[0,2].get_yaxis().set_ticklabels(["", "", ""]),
axs[0,3].get_yaxis().set_ticklabels(["", "", ""])
axs[0,4].get_yaxis().set_ticklabels(["", "", ""])

axs[1,1].get_yaxis().set_ticklabels(["", "", ""])
axs[1,2].get_yaxis().set_ticklabels(["", "", ""])
axs[1,2].get_yaxis().set_ticklabels(["", "", ""])
axs[1,3].get_yaxis().set_ticklabels(["", "", ""])
axs[1,4].get_yaxis().set_ticklabels(["", "", ""])



axs[0,0].set_ylim(0,110)
axs[0,1].set_ylim(0, 110)
axs[0,2].set_ylim(0,110)
axs[0,3].set_ylim(0,110)
axs[0,4].set_ylim(0,110)

axs[1,0].set_ylim(0, 110)
axs[1,1].set_ylim(0,110)
axs[1,2].set_ylim(0,110)
axs[1,2].set_ylim(0,110)
axs[1,3].set_ylim(0,110)
axs[1,4].set_ylim(0,110)




#Place Letter above each bar
#Letters = ("1", "2", "3", "4", "5", "6")


rcParams['font.size']=40
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42


axs[0,0].spines['top'].set_visible(False)
axs[0,1].spines['top'].set_visible(False)
axs[0,2].spines['top'].set_visible(False)
axs[0,3].spines['top'].set_visible(False)
axs[0,4].spines['top'].set_visible(False)

axs[1,0].spines['top'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[1,3].spines['top'].set_visible(False)
axs[1,4].spines['top'].set_visible(False)

axs[0,0].spines['right'].set_visible(False)
axs[0,1].spines['right'].set_visible(False)
axs[0,2].spines['right'].set_visible(False)
axs[0,3].spines['right'].set_visible(False)
axs[0,4].spines['right'].set_visible(False)

axs[1,0].spines['right'].set_visible(False)
axs[1,1].spines['right'].set_visible(False)
axs[1,2].spines['right'].set_visible(False)
axs[1,3].spines['right'].set_visible(False)
axs[1,4].spines['right'].set_visible(False)














# Add some text for labels, title and custom x-axis tick labels, etc.

#axs[1,0].set_title('Indolicidin (3:1)', size=20)
#axs[2,0].set_title('TemporinA (3:1)', size=20)

#axs[1,1].set_title('Indolicidin (9:1)', size=20)
#axs[2,1].set_title('TemporinA (9:1)', size=20)
#axs[0,0].set_ylabel('Relative Intensity', size=12)
#axs[1,0].set_ylabel('Relative Intensity', size=12)
#axs[2,0].set_ylabel('Relative Intensity', size=12)
#axs[0,1].yaxis.set_label_position("right")
#axs[0,1].set_ylabel('Relative Intensity', size=24, rotation=0)
#axs[1,1].yaxis.set_label_position("right")
#axs[1,1].yaxis.set_label_coords(-1.42, 0.5)
#axs[1,1].set_ylabel('Relative Intensity', size=24, rotation=0)
#axs[2,1].yaxis.set_label_position("right")
#axs[2,1].yaxis.set_label_coords(-1.42, 0.5)
#axs[2,1].set_ylabel('Relative Intensity', size=24, rotation=0)

axs[0,0].get_yaxis().set_visible(True)
axs[0,1].get_yaxis().set_visible(True)
axs[0,2].get_yaxis().set_visible(True)
axs[0,3].get_yaxis().set_visible(True)
axs[0,4].get_yaxis().set_visible(True)

axs[1,0].get_yaxis().set_visible(True)
axs[1,1].get_yaxis().set_visible(True)
axs[1,2].get_yaxis().set_visible(True)
axs[1,3].get_yaxis().set_visible(True)
axs[1,4].get_yaxis().set_visible(True)



#axs[1,0].set_yaxis_ticks().set_visible(True)





#axs[5,1].get_xaxis().set_label('Oligomeric State')


#axs[5,1].set_xaxisticks('1', '2', '3', '4', '5', '6')
axs[1,0].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1,1].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 3].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 4].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)

#axsaxs[1, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)[0,0].get_yaxis().set_ticklabels(['0', '100'])


axs[0,0].set_title('0 μM AMT', size=24)
axs[0,1].set_title('19 μM AMT', size=24)
axs[0,2].set_title('37 μM AMT', size=24)
axs[0,3].set_title('75 μM AMT', size=24)
axs[0,4].set_title('150 μM AMT', size=24)


axs[0,0].set_ylabel('WT', size=30, rotation=90)
axs[1,0].set_ylabel('S31N', size=30, rotation=90)

#axs[5,0].set.xlabel = 'Mass (kDa)'

#Set Y-axis range and tick mark spacing

#axs[0,0].set_ylim(0, 100)
#axs[1,0].set_yticks(np.arange(0, 100, step=2))
#axs[1,0].set_ylim(0, 100)
#axs[2,0].set_yticks(np.arange(0, 100, step=2))
#axs[2,0].set_ylim(0, 100)
#axs[0,1].set_yticks(np.arange(0,100, step=2))
#axs[0,1].set_ylim(0, 100)
#axs[1,1].set_yticks(np.arange(0, 100, step=2))
#axs[1,1].set_ylim(0, 100)
#axs[2,1].set_yticks(np.arange(0,100, step=2))
#axs[2,1].set_ylim(0,100)

#Turn off x-axis labels and tick marks
#axs[0,0].get_xaxis().set_visible(False)
#axs[1,0].get_xaxis().set_visible(False)
#axs[2,0].get_xaxis().set_visible(True)
#axs[0,1].get_xaxis().set_visible(False)
#axs[1,1].get_xaxis().set_visible(False)
#axs[2,1].get_xaxis().set_visible(True)
#axs[0,0].spines['top'].set_visible(False)
#axs[1,0].spines['top'].set_visible(False)
#axs[2,0].spines['top'].set_visible(False)
#axs[0,1].spines['top'].set_visible(False)
#axs[1,1].spines['top'].set_visible(False)
#axs[2,1].spines['top'].set_visible(False)
#axs[0,0].spines['right'].set_visible(False)
#axs[1,0].spines['right'].set_visible(False)
#axs[2,0].spines['right'].set_visible(False)
#axs[0,1].spines['right'].set_visible(False)
#axs[1,1].spines['right'].set_visible(False)
#axs[2,1].spines['right'].set_visible(False)

#axs[2,0].tick_params('x', labelsize=20)
#axs[2,1].tick_params('x', labelsize=20)
#axs[0,0].tick_params('y', labelsize=16)
#axs[1,0].tick_params('y', labelsize=16)
#axs[2,0].tick_params('y', labelsize=16)
#axs[0,1].tick_params('y', labelsize=16)
#axs[1,1].tick_params('y', labelsize=16)
#axs[2,1].tick_params('y', labelsize=16)
#axs[2,0].tick_params(axis='x', length=0)
#axs[2,1].tick_params(axis='x', length=0)

plt.savefig("plotter.pdf", transparent="True")
plt.show()
