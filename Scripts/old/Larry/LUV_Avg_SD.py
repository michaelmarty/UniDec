import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import tkinter
from tkinter import filedialog as fd

#files = ["Z:\\Group Share\\Larry Walker\\Fluorescence\\Time course experiments\\DMPG\\Normalized_Data\\3to1_indolicidin_DMPG_Averages.csv",
#         "Z:\\Group Share\\Larry Walker\\Fluorescence\\Time course experiments\\DMPG\\Normalized_Data\\9to1_indolicidin_DMPG_Averages.csv",
#"Z:\\Group Share\\Larry Walker\\Fluorescence\\Time course experiments\\DMPG\\Normalized_Data\\18to1_indolicidin_DMPG_Averages.csv",
#"Z:\\Group Share\\Larry Walker\\Fluorescence\\Time course experiments\\DMPG\\Normalized_Data\\10ul_MeOH_DMPG_Averages.csv"]

root = tkinter.Tk()
files=fd.askopenfilenames(initialdir="D:\Group Share\Fluorescence\Time course experiments", title="Select Files", filetypes=(("csv files", "*.csv"), ("txt files", "*.txt")))
root.destroy()

#colors=["r", "b", "g", "k"]

for i,f in enumerate(files):
    df = pd.read_csv(f)

    stop=-20
    times = df["Time"].to_numpy()[:stop]

    for k in df.keys():
        if "Mean" in k:
            means = df[k].to_numpy()[:stop]
        if "Std" in k:
            stds = df[k].to_numpy()[:stop]/np.sqrt(3)

    #means = df.to_numpy()[:stop,3]
    #stds = means * 0

    #color=colors[i]
    plt.plot(times, means, linestyle="-")
    plt.fill_between(times, means-stds,means+stds, alpha=0.1) #color=color, alpha=0.1)
    plt.legend(files)

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
plt.xlabel("Time (mins)", loc='center')
plt.xlim(0, 30.02)
plt.ylabel("Normalized Intensity", loc='center')
plt.ylim(0, 1.01)
plt.yticks(np.arange(0, 1.01, step=0.25))
#plt.legend#(loc='center left', bbox_to_anchor=(1, 0.5)) # putting legend outside of plot
plt.show()
