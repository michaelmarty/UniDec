import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from matplotlib import rcParams
from matplotlib import cm
import itertools


path = "Z:\Group Share\Julia\M2 protein\wildtype/NEW PURIFICATION/PHGRADIENT"
os.chdir(path)



#Create empty lists for the names of peptides and the averages that will be plotted
ph_4=()
ph_5=()
ph_6=()
ph_7=()
ph_8=()
ph_9=()
Averages_c8e4_4=()
Averages_c8e4_5=()
Averages_c8e4_6=()
Averages_c8e4_7=()
Averages_c8e4_8=()
Averages_c8e4_9=()
Averages_ddm_4=()
Averages_ddm_5=()
Averages_ddm_6=()
Averages_ddm_7=()
Averages_ddm_8=()
Averages_ddm_9=()
Averages_ddm_4=()
Averages_ddm_5=()
Averages_ddm_7=()
Averages_ddm_9=()
Errors_c8e4_4=()
Errors_c8e4_5=()
Errors_c8e4_6=()
Errors_c8e4_7=()
Errors_c8e4_8=()
Errors_c8e4_9=()
Errors_ddm_4=()
Errors_ddm_5=()
Errors_ddm_6=()
Errors_ddm_7=()
Errors_ddm_8=()
Errors_ddm_9=()
Errors_ddm_4=()
Errors_ddm_5=()
Errors_ddm_7=()
Errors_ddm_9=()


#Load files and append values to empty lists created above

ph_4 = np.append(Averages_c8e4_4, np.loadtxt("c8e4_ph4.txt", dtype=str, delimiter='\t'))
ph_5 = np.append(Averages_c8e4_5, np.loadtxt("c8e4_ph5.txt", dtype=str, delimiter='\t'))
ph_6 = np.append(Averages_c8e4_6, np.loadtxt("c8e4_ph6.txt", dtype=str, delimiter='\t'))
ph_7 = np.append(Averages_c8e4_7, np.loadtxt("c8e4_ph7.txt", dtype=str, delimiter='\t'))
ph_8 = np.append(Averages_c8e4_8, np.loadtxt("c8e4_ph8.txt", dtype=str, delimiter='\t'))
ph_9 = np.append(Averages_c8e4_9, np.loadtxt("c8e4_ph9.txt", dtype=str, delimiter='\t'))

Averages_c8e4_4 = np.append(Averages_c8e4_4, np.loadtxt("c8e4_ph4.txt", delimiter='\t'))
Averages_c8e4_5 = np.append(Averages_c8e4_5, np.loadtxt("c8e4_ph5.txt", delimiter='\t'))
Averages_c8e4_6 = np.append(Averages_c8e4_6, np.loadtxt("c8e4_ph6.txt", delimiter='\t'))
Averages_c8e4_7 = np.append(Averages_c8e4_7, np.loadtxt("c8e4_ph7.txt", delimiter='\t'))
Averages_c8e4_8 = np.append(Averages_c8e4_8, np.loadtxt("c8e4_ph8.txt", delimiter='\t'))
Averages_c8e4_9 = np.append(Averages_c8e4_9, np.loadtxt("c8e4_ph9.txt", delimiter='\t'))
Errors_c8e4_4 = np.append(Errors_c8e4_4, np.loadtxt("error_c8e4_ph4.txt", delimiter='\t'))
Errors_c8e4_5 = np.append(Errors_c8e4_5, np.loadtxt("error_c8e4_ph5.txt", delimiter='\t'))
Errors_c8e4_6 = np.append(Errors_c8e4_6, np.loadtxt("error_c8e4_ph6.txt", delimiter='\t'))
Errors_c8e4_7 = np.append(Errors_c8e4_7, np.loadtxt("error_c8e4_ph7.txt", delimiter='\t'))
Errors_c8e4_8 = np.append(Errors_c8e4_8, np.loadtxt("error_c8e4_ph8.txt", delimiter='\t'))
Errors_c8e4_9 = np.append(Errors_c8e4_9, np.loadtxt("error_c8e4_ph9.txt", delimiter='\t'))

Averages_ddm_4 = np.append(Averages_ddm_4, np.loadtxt("ddm_4.txt", delimiter='\t'))
Averages_ddm_5 = np.append(Averages_ddm_5, np.loadtxt("ddm_5.txt", delimiter='\t'))
Averages_ddm_6 = np.append(Averages_ddm_6, np.loadtxt("ddm_6.txt", delimiter='\t'))
Averages_ddm_7 = np.append(Averages_ddm_7, np.loadtxt("ddm_7.txt", delimiter='\t'))
Averages_ddm_8 = np.append(Averages_ddm_8, np.loadtxt("ddm_8.txt", delimiter='\t'))
Averages_ddm_9 = np.append(Averages_ddm_9, np.loadtxt("ddm_9.txt", delimiter='\t'))
Errors_ddm_4 = np.append(Errors_ddm_4, np.loadtxt("error_ddm_4.txt", delimiter='\t'))
Errors_ddm_5 = np.append(Errors_ddm_5, np.loadtxt("error_ddm_5.txt", delimiter='\t'))
Errors_ddm_6 = np.append(Errors_ddm_6, np.loadtxt("error_ddm_6.txt", delimiter='\t'))
Errors_ddm_7 = np.append(Errors_ddm_7, np.loadtxt("error_ddm_7.txt", delimiter='\t'))
Errors_ddm_8 = np.append(Errors_ddm_8, np.loadtxt("error_ddm_8.txt", delimiter='\t'))
Errors_ddm_9 = np.append(Errors_ddm_9, np.loadtxt("error_ddm_9.txt", delimiter='\t'))

Averages_ddm_4 = np.append(Averages_ddm_4, np.loadtxt("ddm_4.txt", delimiter='\t'))
Averages_ddm_5 = np.append(Averages_ddm_5, np.loadtxt("ddm_5.txt", delimiter='\t'))
Averages_ddm_7 = np.append(Averages_ddm_7, np.loadtxt("ddm_7.txt", delimiter='\t'))
Averages_ddm_9 = np.append(Averages_ddm_9, np.loadtxt("ddm_9.txt", delimiter='\t'))
Errors_ddm_4 = np.append(Errors_ddm_4, np.loadtxt("error_ddm_4.txt", delimiter='\t'))
Errors_ddm_5 = np.append(Errors_ddm_5, np.loadtxt("error_ddm_5.txt", delimiter='\t'))
Errors_ddm_7 = np.append(Errors_ddm_7, np.loadtxt("error_ddm_7.txt", delimiter='\t'))
Errors_ddm_9 = np.append(Errors_ddm_9, np.loadtxt("error_ddm_9.txt", delimiter='\t'))




#Set colormap that will be used in bar chart

cmap = matplotlib.cm.get_cmap('viridis_r')
colors = cmap(np.arange(0, len(ph_4)) / (len(ph_4) - 1))

#plt.figure()
#barlist=plt.bar([1,2,3,4,5,6], [1,2,3,4,5,6])
#barlist[0].setcolor('r')
#barlist[1].setcolor('y')
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1),
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1))

#plt.bar = colors = ['black', 'red', 'green', 'blue', 'cyan', 'magenta']




#set the number or subplots in the figure
fig, axs = plt.subplots(3, 6, sharex = 'all')

#Plot the bar charts as subplots
rects1 = axs[0,0].bar(ph_4, Averages_c8e4_4, yerr=Errors_c8e4_4, capsize=3, color=colors)
rects2 = axs[0,1].bar(ph_5, Averages_c8e4_5, yerr=Errors_c8e4_5, capsize=3, color=colors)
rects3 = axs[0,2].bar(ph_6, Averages_c8e4_6, yerr=Errors_c8e4_6, capsize=3, color=colors)
rects4 = axs[0,3].bar(ph_7, Averages_c8e4_7, yerr=Errors_c8e4_7, capsize=3, color=colors)
rects5 = axs[0,4].bar(ph_8, Averages_c8e4_8, yerr=Errors_c8e4_8, capsize=3, color=colors)
rects6 = axs[0,5].bar(ph_9, Averages_c8e4_9, yerr=Errors_c8e4_9, capsize=3, color=colors)
rects7 = axs[1,0].bar(ph_4, Averages_ddm_4, yerr=Errors_ddm_4, capsize=3, color=colors)
rects8 = axs[1,1].bar(ph_5, Averages_ddm_5, yerr=Errors_ddm_5, capsize=3, color=colors)
rects9 = axs[1,2].bar(ph_6, Averages_ddm_6, yerr=Errors_ddm_6, capsize=3, color=colors)
rects10 = axs[1,3].bar(ph_7, Averages_ddm_7, yerr=Errors_ddm_7, capsize=3, color=colors)
rects11 = axs[1,4].bar(ph_8, Averages_ddm_8, yerr=Errors_ddm_8, capsize=3, color=colors)
rects12 = axs[1,5].bar(ph_9, Averages_ddm_9, yerr=Errors_ddm_9, capsize=3, color=colors)
rects13 = axs[2,0].bar(ph_4, Averages_ddm_4, yerr=Errors_ddm_4, capsize=3, color=colors)
rects14 = axs[2,1].bar(ph_5, Averages_ddm_5, yerr=Errors_ddm_5, capsize=3, color=colors)
rects15 = axs[2,3].bar(ph_7, Averages_ddm_7, yerr=Errors_ddm_7, capsize=3, color=colors)
rects16 = axs[2,5].bar(ph_9, Averages_ddm_9, yerr=Errors_ddm_9, capsize=3, color=colors)




axs[0, 0].get_yaxis().set_ticks([0, 50, 100])
axs[0, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[1, 0].get_yaxis().set_ticks([0, 50, 100])
axs[1, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[2, 0].get_yaxis().set_ticks([0, 50, 100])
axs[2, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),

axs[0, 1].get_yaxis().set_ticks([0, 50, 100])
axs[0, 2].get_yaxis().set_ticks([0, 50, 100])
axs[0, 3].get_yaxis().set_ticks([0, 50, 100])
axs[0, 4].get_yaxis().set_ticks([0, 50, 100])
axs[0, 5].get_yaxis().set_ticks([0, 50, 100])
axs[1, 1].get_yaxis().set_ticks([0, 50, 100])
axs[1, 2].get_yaxis().set_ticks([0, 50, 100])
axs[1, 3].get_yaxis().set_ticks([0, 50, 100])
axs[1, 4].get_yaxis().set_ticks([0, 50, 100])
axs[1, 5].get_yaxis().set_ticks([0, 50, 100])
axs[1, 2].get_yaxis().set_ticks([0, 50, 100])

axs[0, 1].get_yaxis().set_ticklabels(["", '', ""]),
axs[0, 2].get_yaxis().set_ticklabels(["", '', ""]),
axs[0, 3].get_yaxis().set_ticklabels(["", '', ""]),
axs[0, 4].get_yaxis().set_ticklabels(["", '', ""]),
axs[0, 5].get_yaxis().set_ticklabels(["", '', ""]),
axs[1, 1].get_yaxis().set_ticklabels(["", '', ""]),
axs[1, 2].get_yaxis().set_ticklabels(["", '', ""]),
axs[1, 3].get_yaxis().set_ticklabels(["", '', ""]),
axs[1, 4].get_yaxis().set_ticklabels(["", '', ""]),
axs[1, 5].get_yaxis().set_ticklabels(["", '', ""]),
axs[2, 1].get_yaxis().set_ticklabels(["", '', ""]),
axs[2, 2].get_yaxis().set_ticklabels(["", '', ""]),
axs[2, 3].get_yaxis().set_ticklabels(["", '', ""]),
axs[2, 4].get_yaxis().set_ticklabels(["", '', ""]),
axs[2, 5].get_yaxis().set_ticklabels(["", '', ""]),






axs[0,0].set_ylim(0,115)
axs[1,0].set_ylim(0, 115)
axs[2,0].set_ylim(0,115)
axs[0,1].set_ylim(0,115)
axs[0,2].set_ylim(0,115)
axs[0,3].set_ylim(0,115)
axs[0,4].set_ylim(0, 115)
axs[0,5].set_ylim(0,115)
axs[1,1].set_ylim(0,115)
axs[1,2].set_ylim(0,115)
axs[1,3].set_ylim(0,115)
axs[1,4].set_ylim(0,115)
axs[1,5].set_ylim(0,115)
axs[2,1].set_ylim(0,115)
axs[2,2].set_ylim(0,115)
axs[2,3].set_ylim(0,115)
axs[2,4].set_ylim(0,115)
axs[2,5].set_ylim(0,115)



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
axs[0,5].spines['top'].set_visible(False)
axs[1,0].spines['top'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[1,3].spines['top'].set_visible(False)
axs[1,4].spines['top'].set_visible(False)
axs[1,5].spines['top'].set_visible(False)
axs[2,0].spines['top'].set_visible(False)
axs[2,1].spines['top'].set_visible(False)
axs[2,2].spines['top'].set_visible(False)
axs[2,3].spines['top'].set_visible(False)
axs[2,4].spines['top'].set_visible(False)
axs[2,5].spines['top'].set_visible(False)


axs[0,0].spines['right'].set_visible(False)
axs[0,1].spines['right'].set_visible(False)
axs[0,2].spines['right'].set_visible(False)
axs[0,3].spines['right'].set_visible(False)
axs[0,4].spines['right'].set_visible(False)
axs[0,5].spines['right'].set_visible(False)
axs[1,0].spines['right'].set_visible(False)
axs[1,1].spines['right'].set_visible(False)
axs[1,2].spines['right'].set_visible(False)
axs[1,3].spines['right'].set_visible(False)
axs[1,4].spines['right'].set_visible(False)
axs[1,5].spines['right'].set_visible(False)
axs[2,0].spines['right'].set_visible(False)
axs[2,1].spines['right'].set_visible(False)
axs[2,2].spines['right'].set_visible(False)
axs[2,3].spines['right'].set_visible(False)
axs[2,4].spines['right'].set_visible(False)
axs[2,5].spines['right'].set_visible(False)



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
axs[0,5].get_yaxis().set_visible(True)
axs[1,0].get_yaxis().set_visible(True)
axs[1,1].get_yaxis().set_visible(True)
axs[1,2].get_yaxis().set_visible(True)
axs[1,3].get_yaxis().set_visible(True)
axs[1,4].get_yaxis().set_visible(True)
axs[1,5].get_yaxis().set_visible(True)
axs[2,0].get_yaxis().set_visible(True)
axs[2,1].get_yaxis().set_visible(True)
axs[2,2].get_yaxis().set_visible(True)
axs[2,3].get_yaxis().set_visible(True)
axs[2,4].get_yaxis().set_visible(True)
axs[2,5].get_yaxis().set_visible(True)

#axs[1,0].set_yaxis_ticks().set_visible(True)



axs[0,0].get_yaxis().set_label('Relative Intensity')





#axs[5,1].set_xaxisticks('1', '2', '3', '4', '5', '6')
axs[2, 1].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[2, 0].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[2, 3].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[2, 4].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[2, 5].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[2, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)

#axs[0,0].get_yaxis().set_ticklabels(['0', '100'])


axs[0,0].set_title('pH 4', size=30)
axs[0,1].set_title('pH 5', size=30)
axs[0,2].set_title('pH 6', size=30)
axs[0,3].set_title('pH 7', size=30)
axs[0,4].set_title('pH 8', size=30)
axs[0,5].set_title('pH 9', size=30)


axs[0,0].set_ylabel('C8E4', size=24, rotation=90)
axs[1,0].set_ylabel('ddm', size=24, rotation=90)
axs[2,0].set_ylabel('DDM', size=24, rotation=90)


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