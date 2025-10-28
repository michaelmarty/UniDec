import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as stats
from matplotlib import rcParams
from matplotlib import cm
import itertools


path = "Z:\\Group Share\\Julia\\M2 protein\\wildtype\\NEW PURIFICATION\\DETERGENT SCREEN PEAK INTEGRATION"
os.chdir(path)



#Create empty lists for the names of peptides and the averages that will be plotted
ph_5=()
ph_7=()
ph_9=()
Averages_c8e4_5=()
Averages_c8e4_7=()
Averages_c8e4_9=()
Averages_ddm_5=()
Averages_ddm_7=()
Averages_ddm_9=()
Averages_og_5=()
Averages_og_7=()
Averages_og_9=()
Averages_ddm_5=()
Averages_ddm_7=()
Averages_ddm_9=()
Averages_dpc_5=()
Averages_dpc_7=()
Averages_dpc_9=()
Averages_lmng_5=()
Averages_lmng_7=()
Averages_lmng_9=()
Errors_c8e4_5=()
Errors_c8e4_7=()
Errors_c8e4_9=()
Errors_ddm_5=()
Errors_ddm_7=()
Errors_ddm_9=()
Errors_og_5=()
Errors_og_7=()
Errors_og_9=()
Errors_ddm_5=()
Errors_ddm_7=()
Errors_ddm_9=()
Errors_dpc_5=()
Errors_dpc_7=()
Errors_dpc_9=()
Errors_lmng_5=()
Errors_lmng_7=()
Errors_lmng_9=()

#Load files and append values to empty lists created above
ph_5 = np.append(Averages_c8e4_5, np.loadtxt("c8e4_ph5.txt", dtype=str, delimiter='\t'))
Averages_c8e4_5 = np.append(Averages_c8e4_5, np.loadtxt("c8e4_ph5.txt", delimiter='\t'))
Errors_c8e4_5 = np.append(Errors_c8e4_5, np.loadtxt("error_c8e4_ph5.txt", delimiter='\t'))
ph_7= np.append(Averages_c8e4_7, np.loadtxt("c8e4_ph7.txt", dtype=str, delimiter='\t'))
Averages_c8e4_7 = np.append(Averages_c8e4_7, np.loadtxt("c8e4_ph7.txt", delimiter='\t'))
Errors_c8e4_7 = np.append(Errors_c8e4_7, np.loadtxt("error_c8e4_ph7.txt", delimiter='\t'))
ph_9 = np.append(Averages_c8e4_9, np.loadtxt("c8e4_ph9.txt", dtype=str, delimiter='\t'))
Averages_c8e4_9 = np.append(Averages_c8e4_9, np.loadtxt("c8e4_ph9.txt", delimiter='\t'))
Errors_c8e4_9 = np.append(Errors_c8e4_9, np.loadtxt("error_c8e4_ph9.txt", delimiter='\t'))
Averages_ddm_5 = np.append(Averages_ddm_5, np.loadtxt("ddm_ph5.txt", delimiter='\t'))
Errors_ddm_5 = np.append(Errors_ddm_5, np.loadtxt("error_ddm_ph5.txt", delimiter='\t'))
Averages_ddm_7 = np.append(Averages_ddm_7, np.loadtxt("ddm_ph7.txt", delimiter='\t'))
Errors_ddm_7 = np.append(Errors_ddm_7, np.loadtxt("error_ddm_ph7.txt", delimiter='\t'))
Averages_ddm_9 = np.append(Averages_ddm_9, np.loadtxt("ddm_ph9.txt", delimiter='\t'))
Errors_ddm_9 = np.append(Errors_ddm_9, np.loadtxt("error_ddm_ph9.txt", delimiter='\t'))
Averages_og_5 = np.append(Averages_og_5, np.loadtxt("og_ph5.txt", delimiter='\t'))
Errors_og_5 = np.append(Errors_og_5, np.loadtxt("error_og_ph5.txt", delimiter='\t'))
Averages_og_7 = np.append(Averages_og_7, np.loadtxt("og_ph7.txt", delimiter='\t'))
Errors_og_7 = np.append(Errors_og_7, np.loadtxt("error_og_ph7.txt", delimiter='\t'))
Averages_og_9 = np.append(Averages_og_9, np.loadtxt("og_ph9.txt", delimiter='\t'))
Errors_og_9 = np.append(Errors_og_9, np.loadtxt("error_og_ph9.txt", delimiter='\t'))
Averages_ddm_5 = np.append(Averages_ddm_5, np.loadtxt("ddm_ph5.txt", delimiter='\t'))
Errors_ddm_5 = np.append(Errors_ddm_5, np.loadtxt("error_ddm_ph5.txt", delimiter='\t'))
Averages_ddm_7 = np.append(Averages_ddm_7, np.loadtxt("ddm_ph7.txt", delimiter='\t'))
Errors_ddm_7 = np.append(Errors_ddm_7, np.loadtxt("error_ddm_ph7.txt", delimiter='\t'))
Averages_ddm_9 = np.append(Averages_ddm_9, np.loadtxt("ddm_ph9.txt", delimiter='\t'))
Errors_ddm_9 = np.append(Errors_ddm_9, np.loadtxt("error_ddm_ph9.txt", delimiter='\t'))
Averages_dpc_5 = np.append(Averages_dpc_5, np.loadtxt("dpc_ph5.txt", delimiter='\t'))
Errors_dpc_5 = np.append(Errors_dpc_5, np.loadtxt("error_dpc_ph5.txt", delimiter='\t'))
Averages_dpc_7 = np.append(Averages_dpc_7, np.loadtxt("dpc_ph7.txt", delimiter='\t'))
Errors_dpc_7 = np.append(Errors_dpc_7, np.loadtxt("error_dpc_ph7.txt", delimiter='\t'))
Averages_dpc_9 = np.append(Averages_dpc_9, np.loadtxt("dpc_ph9.txt", delimiter='\t'))
Errors_dpc_9 = np.append(Errors_dpc_9, np.loadtxt("error_dpc_ph9.txt", delimiter='\t'))
Averages_lmng_5 = np.append(Averages_lmng_5, np.loadtxt("lmng_ph5.txt", delimiter='\t'))
Errors_lmng_5 = np.append(Errors_lmng_5, np.loadtxt("error_lmng_ph5.txt", delimiter='\t'))
Averages_lmng_7 = np.append(Averages_lmng_5, np.loadtxt("lmng_ph7.txt", delimiter='\t'))
Errors_lmng_7 = np.append(Errors_lmng_7, np.loadtxt("error_lmng_ph7.txt", delimiter='\t'))
Averages_lmng_9 = np.append(Averages_lmng_9, np.loadtxt("lmng_ph9.txt", delimiter='\t'))
Errors_lmng_9 = np.append(Errors_lmng_9, np.loadtxt("error_lmng_ph9.txt", delimiter='\t'))



#Set colormap that will be used in bar chart

cmap = matplotlib.cm.get_cmap('viridis_r')
colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1))

#plt.figure()
#barlist=plt.bar([1,2,3,4,5,6], [1,2,3,4,5,6])
#barlist[0].setcolor('r')
#barlist[1].setcolor('y')
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1),
#colors = cmap(np.arange(0, len(ph_5)) / (len(ph_5) - 1))

#plt.bar = colors = ['black', 'red', 'green', 'blue', 'cyan', 'magenta']




#set the number or subplots in the figure
fig, axs = plt.subplots(6, 3, sharex = 'all')

#Plot the bar charts as subplots
rects1 = axs[0,0].bar(ph_5, Averages_c8e4_5, yerr=Errors_c8e4_5, capsize=3, color=colors)
rects2 = axs[1,0].bar(ph_5, Averages_ddm_5, yerr=Errors_ddm_5, capsize=3, color=colors)
rects3 = axs[2,0].bar(ph_5, Averages_og_5, yerr=Errors_og_5,  capsize=3, color=colors)
rects4 = axs[3,0].bar(ph_5, Averages_ddm_5, yerr=Errors_ddm_5, capsize=3, color=colors)
rects5 = axs[4,0].bar(ph_5, Averages_dpc_5, yerr=Errors_dpc_5, capsize=3, color=colors)
rects6 = axs[5,0].bar(ph_5, Averages_lmng_5, yerr=Errors_lmng_5,  capsize=3, color=colors)
rects7 = axs[0,1].bar(ph_7, Averages_c8e4_7, yerr=Errors_c8e4_7, capsize=3, color=colors)
rects8 = axs[1,1].bar(ph_7, Averages_ddm_7, yerr=Errors_ddm_7, capsize=3, color=colors)
rects9 = axs[2,1].bar(ph_7, Averages_og_7, yerr=Errors_og_7, capsize=3, color=colors)
rects10 = axs[3,1].bar(ph_7, Averages_ddm_7, yerr=Errors_ddm_7, capsize=3, color=colors)
rects11 = axs[4,1].bar(ph_7, Averages_dpc_7, yerr=Errors_dpc_7, capsize=3, color=colors)
#rects12 = axs[5,1].bar(ph_7, Averages_lmng_7, yerr=Errors_lmng_7, capsize=3)
rects13 = axs[0,2].bar(ph_9, Averages_c8e4_9, yerr=Errors_c8e4_9, capsize=3, color=colors)
rects14 = axs[1,2].bar(ph_9, Averages_ddm_9, yerr=Errors_ddm_9, capsize=3, color=colors)
rects15 = axs[2,2].bar(ph_9, Averages_og_9, yerr=Errors_og_9, capsize=3, color=colors)
rects16 = axs[3,2].bar(ph_9, Averages_ddm_9, yerr=Errors_ddm_9, capsize=3, color=colors)
rects17 = axs[4,2].bar(ph_9, Averages_dpc_9, yerr=Errors_dpc_9, capsize=3, color=colors)
rects18 = axs[5,2].bar(ph_9, Averages_lmng_9, yerr=Errors_lmng_9, capsize=3, color=colors)


axs[0, 0].get_yaxis().set_ticks([0, 50, 100])
axs[0, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20)
axs[1, 0].get_yaxis().set_ticks([0, 50, 100])
axs[1, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[2, 0].get_yaxis().set_ticks([0, 50, 100])
axs[2, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[3, 0].get_yaxis().set_ticks([0, 50, 100])
axs[3, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[4, 0].get_yaxis().set_ticks([0, 50, 100])
axs[4, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[5, 0].get_yaxis().set_ticks([0, 50, 100])
axs[5, 0].get_yaxis().set_ticklabels(["0", '%', "100"], size=20),
axs[0, 1].get_yaxis().set_ticks([0, 50, 100])
axs[1, 1].get_yaxis().set_ticks([0, 50, 100])
axs[2, 1].get_yaxis().set_ticks([0, 50, 100])
axs[3, 1].get_yaxis().set_ticks([0, 50, 100])
axs[4, 1].get_yaxis().set_ticks([0, 50, 100])
axs[5, 1].get_yaxis().set_ticks([0, 50, 100])
axs[0, 2].get_yaxis().set_ticks([0, 50, 100])
axs[1, 2].get_yaxis().set_ticks([0, 50, 100])
axs[2, 2].get_yaxis().set_ticks([0, 50, 100])
axs[3, 2].get_yaxis().set_ticks([0, 50, 100])
axs[4, 2].get_yaxis().set_ticks([0, 50, 100])
axs[5, 2].get_yaxis().set_ticks([0, 50, 100])
axs[5, 1].get_yaxis().set_ticklabels(["", '', ""]),
axs[0, 1].get_yaxis().set_ticklabels(["", '', ""]),
axs[1,1].get_yaxis().set_ticklabels(["", "", ""]),
axs[2,1].get_yaxis().set_ticklabels(["", "", ""])
axs[3,1].get_yaxis().set_ticklabels(["", "", ""])
axs[4,1].get_yaxis().set_ticklabels(["", "", ""])
axs[5,1].get_yaxis().set_ticklabels(["", "", ""])
axs[0,2].get_yaxis().set_ticklabels(["", "", ""])
axs[1,2].get_yaxis().set_ticklabels(["", "", ""])
axs[2,2].get_yaxis().set_ticklabels(["", "", ""])
axs[3,2].get_yaxis().set_ticklabels(["", "", ""])
axs[4,2].get_yaxis().set_ticklabels(["", "", ""])
axs[5,2].get_yaxis().set_ticklabels(["", "", ""])


axs[0,0].set_ylim(0,110)
axs[1,0].set_ylim(0, 110)
axs[2,0].set_ylim(0,110)
axs[3,0].set_ylim(0,110)
axs[4,0].set_ylim(0,110)
axs[5,0].set_ylim(0,110)
axs[0,1].set_ylim(0, 110)
axs[1,1].set_ylim(0,110)
axs[2,1].set_ylim(0,110)
axs[3,1].set_ylim(0,110)
axs[4,1].set_ylim(0,110)
axs[5,1].set_ylim(0,110)
axs[0,2].set_ylim(0,110)
axs[1,2].set_ylim(0,110)
axs[2,2].set_ylim(0,110)
axs[3,2].set_ylim(0,110)
axs[4,2].set_ylim(0,110)
axs[5,2].set_ylim(0,110)



#Place Letter above each bar
#Letters = ("1", "2", "3", "4", "5", "6")


rcParams['font.size']=40
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

plt.title('C8E4')

axs[0,0].spines['top'].set_visible(False)
axs[0,2].spines['top'].set_visible(False)
axs[1,0].spines['right'].set_visible(False)
axs[3,0].spines['right'].set_visible(False)
axs[2,0].spines['right'].set_visible(False)
axs[4,0].spines['right'].set_visible(False)
axs[5,0].spines['right'].set_visible(False)
axs[0,1].spines['right'].set_visible(False)
axs[1,1].spines['right'].set_visible(False)
axs[2,1].spines['right'].set_visible(False)
axs[3,1].spines['right'].set_visible(False)
axs[4,1].spines['right'].set_visible(False)
axs[5,1].spines['right'].set_visible(False)
axs[2,2].spines['top'].set_visible(False)
axs[5,2].spines['top'].set_visible(False)
axs[1,0].spines['top'].set_visible(False)
axs[2,0].spines['top'].set_visible(False)
axs[3,0].spines['top'].set_visible(False)
axs[4,0].spines['top'].set_visible(False)
axs[5,0].spines['top'].set_visible(False)
axs[0,1].spines['top'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[2,1].spines['top'].set_visible(False)
axs[3,1].spines['top'].set_visible(False)
axs[4,1].spines['top'].set_visible(False)
axs[5,1].spines['top'].set_visible(False)
axs[0,2].spines['right'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[3,2].spines['top'].set_visible(False)
axs[4,2].spines['top'].set_visible(False)
axs[5,1].spines['top'].set_visible(False)
axs[0,0].spines['right'].set_visible(False)
axs[1,0].spines['right'].set_visible(False)
axs[3,0].spines['right'].set_visible(False)
axs[2,0].spines['right'].set_visible(False)
axs[4,0].spines['right'].set_visible(False)
axs[5,0].spines['right'].set_visible(False)
axs[0,1].spines['right'].set_visible(False)
axs[1,1].spines['right'].set_visible(False)
axs[2,1].spines['right'].set_visible(False)
axs[3,1].spines['right'].set_visible(False)
axs[4,1].spines['right'].set_visible(False)
axs[5,1].spines['right'].set_visible(False)
axs[0,0].spines['right'].set_visible(False)
axs[1,0].spines['right'].set_visible(False)
axs[3,0].spines['right'].set_visible(False)
axs[2,0].spines['right'].set_visible(False)
axs[4,0].spines['right'].set_visible(False)
axs[5,0].spines['right'].set_visible(False)
axs[0,2].spines['right'].set_visible(False)
axs[1,2].spines['right'].set_visible(False)
axs[2,2].spines['right'].set_visible(False)
axs[3,2].spines['right'].set_visible(False)
axs[4,2].spines['right'].set_visible(False)
axs[5,2].spines['right'].set_visible(False)



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
axs[1,0].get_yaxis().set_visible(True)
axs[2,0].get_yaxis().set_visible(True)
axs[3,0].get_yaxis().set_visible(True)
axs[4,0].get_yaxis().set_visible(True)
axs[5,0].get_yaxis().set_visible(True)
axs[0,1].get_yaxis().set_visible(True)
axs[1,1].get_yaxis().set_visible(True)
axs[2,1].get_yaxis().set_visible(True)
axs[3,1].get_yaxis().set_visible(True)
axs[4,1].get_yaxis().set_visible(True)
axs[5,1].get_yaxis().set_visible(True)
axs[0,2].get_yaxis().set_visible(True)
axs[1,2].get_yaxis().set_visible(True)
axs[2,2].get_yaxis().set_visible(True)
axs[3,2].get_yaxis().set_visible(True)
axs[4,2].get_yaxis().set_visible(True)
axs[5,2].get_yaxis().set_visible(True)

#axs[1,0].set_yaxis_ticks().set_visible(True)



axs[0,0].get_yaxis().set_label('Relative Intensity')


axs[5,1].get_xaxis().set_label('Oligomeric State')


#axs[5,1].set_xaxisticks('1', '2', '3', '4', '5', '6')
axs[5, 1].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[5, 0].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
axs[5, 2].get_xaxis().set_ticklabels(["1", '2', '3', '4', '5', "6"], size=20)
#axs[0,0].get_yaxis().set_ticklabels(['0', '100'])


axs[0,0].set_title('pH 5', size=24)
axs[0,1].set_title('pH 7', size=24)
axs[0,2].set_title('pH 9', size=24)


axs[0,0].set_ylabel('C8E4', size=20, rotation=90)
axs[1,0].set_ylabel('ddm', size=20, rotation=90)
axs[2,0].set_ylabel('OG', size=20, rotation=90)
axs[3,0].set_ylabel('DDM', size=20, rotation=90)
axs[4,0].set_ylabel('DPC', size=20, rotation=90)
axs[5,0].set_ylabel('LMNG', size=20, rotation=90)

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