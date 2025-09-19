import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from matplotlib import rcParams
from matplotlib import cm
import itertools


path = "Z:\Julia\M2 protein\wildtype/NEW PURIFICATION/210121/ddm - text files"
os.chdir(path)



#Create empty lists for the names of peptides and the averages that will be plotted
thirteen=()
twofive=()
fivezero=()
hundred=()
twohundred=()


Averages_ddm_5_13=()
Averages_ddm_5_25=()
Averages_ddm_5_50=()
Averages_ddm_5_100=()
Averages_ddm_5_200=()
Averages_ddm_9_13=()
Averages_ddm_9_25=()
Averages_ddm_9_50=()
Averages_ddm_9_100=()
Averages_ddm_9_200=()

Error_ddm_5_13=()
Error_ddm_5_25=()
Error_ddm_5_50=()
Error_ddm_5_100=()
Error_ddm_5_200=()
Error_ddm_9_13=()
Error_ddm_9_25=()
Error_ddm_9_50=()
Error_ddm_9_100=()
Error_ddm_9_200=()

#Load files and append values to empty lists created above


thirteen = np.append(Averages_ddm_5_13, np.loadtxt("13_ddm_5.txt", dtype=str, delimiter='\t'))
Averages_ddm_5_13 = np.append(Averages_ddm_5_13, np.loadtxt("13_ddm_5.txt", delimiter='\t'))
Error_ddm_5_13= np.append(Error_ddm_5_13, np.loadtxt("error_13_ddm_5.txt", delimiter='\t'))
twofive = np.append(Averages_ddm_5_25, np.loadtxt("25_ddm_5.txt", dtype=str, delimiter='\t'))
Averages_ddm_5_25 = np.append(Averages_ddm_5_25, np.loadtxt("25_ddm_5.txt", delimiter='\t'))
Error_ddm_5_25 = np.append(Error_ddm_5_25, np.loadtxt("error_25_ddm_5.txt",  delimiter='\t'))
fivezero = np.append(Averages_ddm_5_50, np.loadtxt("50_ddm_5.txt", dtype=str, delimiter='\t'))
Averages_ddm_5_50 = np.append(Averages_ddm_5_50, np.loadtxt("50_ddm_5.txt",  delimiter='\t'))
Error_ddm_5_50 = np.append(Error_ddm_5_50, np.loadtxt("error_50_ddm_5.txt", delimiter='\t'))
hundred = np.append(Averages_ddm_5_100, np.loadtxt("100_ddm_5.txt", dtype=str, delimiter='\t'))
Averages_ddm_5_100 = np.append(Averages_ddm_5_100, np.loadtxt("100_ddm_5.txt", delimiter='\t'))
Error_ddm_5_100 = np.append(Error_ddm_5_100, np.loadtxt("error_100_ddm_5.txt", delimiter='\t'))
twohundred = np.append(Averages_ddm_5_200, np.loadtxt("200_ddm_5.txt", dtype=str, delimiter='\t'))
Averages_ddm_5_200 = np.append(Averages_ddm_5_200, np.loadtxt("200_ddm_5.txt", delimiter='\t'))
Error_ddm_5_200 = np.append(Error_ddm_5_200, np.loadtxt("error_200_ddm_5.txt",  delimiter='\t'))

Averages_ddm_9_13 = np.append(Averages_ddm_9_13, np.loadtxt("13_ddm_9.txt", delimiter='\t'))
Error_ddm_9_13 = np.append(Error_ddm_9_13, np.loadtxt("error_13_ddm_9.txt", delimiter='\t'))
Averages_ddm_9_25 = np.append(Averages_ddm_9_25, np.loadtxt("25_ddm_9.txt", delimiter='\t'))
Error_ddm_9_25 = np.append(Error_ddm_9_25, np.loadtxt("error_25_ddm_9.txt", delimiter='\t'))
Averages_ddm_9_50 = np.append(Averages_ddm_9_50, np.loadtxt("50_ddm_9.txt", delimiter='\t'))
Error_ddm_9_50 = np.append(Error_ddm_9_50, np.loadtxt("Error_50_ddm_9.txt", delimiter='\t'))
Averages_ddm_9_100 = np.append(Averages_ddm_9_100, np.loadtxt("100_ddm_9.txt", delimiter='\t'))
Error_ddm_9_100 = np.append(Error_ddm_9_100, np.loadtxt("error_100_ddm_9.txt", delimiter='\t'))
Averages_ddm_9_200 = np.append(Averages_ddm_9_200, np.loadtxt("200_ddm_9.txt", delimiter='\t'))
Error_ddm_9_200 = np.append(Error_ddm_9_200, np.loadtxt("error_200_ddm_9.txt", delimiter='\t'))



#Set colormap that will be used in bar chart

cmap = matplotlib.cm.get_cmap('viridis_r')
colors = cmap(np.arange(0, len(thirteen)) / (len(thirteen) - 1))

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
rects1 = axs[0,0].bar(thirteen, Averages_ddm_5_13, yerr=Error_ddm_5_13, capsize=3, color=colors)
rects2 = axs[0,1].bar(twofive, Averages_ddm_5_25, yerr=Error_ddm_5_25, capsize=3, color=colors)
rects3 = axs[0,2].bar(fivezero, Averages_ddm_5_50, yerr=Error_ddm_5_25, capsize=3, color=colors)
rects4 = axs[0,3].bar(hundred, Averages_ddm_5_100, yerr=Error_ddm_5_100, capsize=3, color=colors)
rects5 = axs[0,4].bar(twohundred, Averages_ddm_5_200, yerr=Error_ddm_5_200, capsize=3, color=colors)
rects6 = axs[1,0].bar(thirteen, Averages_ddm_9_13, yerr=Error_ddm_9_13, capsize=3, color=colors)
rects7 = axs[1,1].bar(twofive, Averages_ddm_9_25, yerr=Error_ddm_9_25, capsize=3, color=colors)
rects8 = axs[1,2].bar(fivezero, Averages_ddm_9_50, yerr=Error_ddm_9_50, capsize=3, color=colors)
rects9 = axs[1,3].bar(hundred, Averages_ddm_9_100, yerr=Error_ddm_9_100, capsize=3, color=colors)
rects10 = axs[1,4].bar(twohundred, Averages_ddm_9_200, yerr=Error_ddm_9_200, capsize=3, color=colors)



fig.suptitle('ddm', fontsize=36)



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



axs[0,0].get_yaxis().set_label('Relative Intensity')


#axs[5,1].get_xaxis().set_label('Oligomeric State')


#axs[5,1].set_xaxisticks('1', '2', '3', '4', '5', '6')
axs[1,0].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1,1].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 3].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[1, 4].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)

#axsaxs[1, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)[0,0].get_yaxis().set_ticklabels(['0', '100'])


axs[0,0].set_title('13 μM', size=24)
axs[0,1].set_title('25 μM', size=24)
axs[0,2].set_title('50 μM', size=24)
axs[0,3].set_title('100 μM', size=24)
axs[0,4].set_title('200 μM', size=24)


axs[0,0].set_ylabel('pH 5', size=30, rotation=90)
axs[1,0].set_ylabel('pH 9', size=30, rotation=90)

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
