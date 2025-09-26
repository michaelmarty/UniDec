import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter import filedialog as fd
import pandas as pd
from matplotlib import rcParams
import scipy.ndimage as nd
import re

#open file dialog box to open multiple files
root = tkinter.Tk()
files=fd.askopenfilenames(initialdir="E:\Group Share\Fluorescence\\Time course experiments\\", title="Select Files", filetypes=(("txt files", "*.txt"),("csv files", "*.csv")))
#split into list of files
files = root.tk.splitlist(files)
root.destroy()
#print(files)

xmax=()
filename=()
for file in files:
    with open(file) as f:
        n = len(list(f)) #get length of each file

    filename = os.path.basename(file)
    filename = os.path.splitext(filename)
    filename = str(filename[0]).split('_')
    filename = filename
    print(filename)
    xdata = np.loadtxt(fname=file, skiprows=4, max_rows=n-5, dtype=str, delimiter='\t', usecols=(0))    #Load file data from first column of file into xdata skipping first 4 rows and last 5 rows
    ydata = np.loadtxt(fname=file, skiprows=4, max_rows=n-5, dtype=str, delimiter='\t', usecols=(1))    #Load file data from second column of file into ydata skipping first 4 rows and last 5 rows
    xdata = pd.to_numeric(xdata)    #convert xdata to numeric values
    xdata[:] = [x / 60 for x in xdata]  #convert xdata from seconds to minutes

    ydata = pd.to_numeric(ydata)    #convert ydata to numeric values
    xmax=np.append(xmax, max(xdata))    #find max xvalue in each file and append to xmax

    for i in ydata:
        ydata = (ydata - ydata.min()) / (ydata.max() - ydata.min()) #normalize intensities
    fig = plt.plot(xdata, ydata, label=filename)    #plot x,y scatter plot with lipid name

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
print(xmax)
plt.xlabel("Time (mins)", loc='center')
plt.xlim(0, max(xmax))
plt.ylabel("Normalized Intensity", loc='center')
plt.ylim(0, 1.01)
plt.yticks(np.arange(0, 1.01, step=0.25))
plt.legend()
plt.show()