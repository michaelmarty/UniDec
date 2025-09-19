import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter import filedialog as fd
import pandas as pd
from matplotlib import rcParams
import csv
import scipy.ndimage as nd
import itertools
import scipy.stats as stats
import re

#open file dialog box to open multiple files
root = tkinter.Tk()
files=fd.askopenfilenames(initialdir="E:\Group Share\Fluorescence\Time course experiments", title="Select Files", filetypes=(("txt files", "*.txt"),("csv files", "*.csv")))
filez=files[0]
directory = os.path.dirname(filez)  #Get directory of first selected file
#Creating a new directory to save normalized data
New_folder = os.path.join(directory, "Normalized_Data")
try:
    os.mkdir(New_folder)
except OSError as error:
    print()
files = root.tk.splitlist(files)    #split into list of files
root.destroy()

newpath=()
xmax=()
filename=()
pep_name=()
pep_ratio=()
lipid_name=()
xvals=()
yvals=()
data_files=()
for file in files:
    with open(file) as f:
        n = len(list(f)) #get length of each file

    filename = os.path.basename(file) #getting the file name
    filename = os.path.splitext(filename)   #splitting the filename
    filename = str(filename[0]).split('_')  #splitting the filename by underscore
    filename = filename
    pep_name = filename[3]
    pep_ratio = filename[2]
    lipid_name = filename[1]
    #data_sets = []
#    print(filename)
    xdata = np.loadtxt(fname=file, skiprows=4, max_rows=n-5, dtype=str, delimiter='\t', usecols=(0))    #Load file data from first column of file into xdata skipping first 4 rows and last 5 rows
    ydata = np.loadtxt(fname=file, skiprows=4, max_rows=n-5, dtype=str, delimiter='\t', usecols=(1))    #Load file data from second column of file into ydata skipping first 4 rows and last 5 rows
    xdata = pd.to_numeric(xdata)    #convert xdata to numeric values
    xdata[:] = [x / 60 for x in xdata]  #convert xdata from seconds to minutes

    ydata = pd.to_numeric(ydata)    #convert ydata to numeric values
    xmax=np.append(xmax, max(xdata))    #find max xvalue in each file and append to xmax
    #yaverage=()
    for i in ydata:
        ydata = (ydata - ydata.min()) / (ydata.max() - ydata.min()) #normalize intensities

    data = np.column_stack([xdata, ydata])
    headers = ['Time', 'Intensity']
    newfile = 'Normalized_' + filename[2] + '_' + filename[3] + '_' + filename[1] + '_' + filename[0] + '.txt'   #creating a new file using parts of original file name
    filepath = os.path.join(New_folder, newfile).replace("/","\\")  #Forcing all of the slashes in directory to be same direction
    newfile = open(filepath, "w+")  #opening the new file for writing
    np.savetxt(filepath, data, header='Time Intensity', comments='')     #saving normalized ydata to new file
    #data_sets.append(pep_ratio + '_' + pep_name)

    New_folder = os.path.join(New_folder).replace("/", "\\")
    data_files = glob.glob(os.path.join(New_folder, 'Normalized_' + pep_ratio + '_' + pep_name + '_' + lipid_name + '*.txt'))    #opening all the normalized data files
    df_final = pd.DataFrame(columns=['Time'])   #creating a pandas data frame to store the results
    #print(data_files)

    #opening and joining the files
    for file in data_files:
        file_df = pd.read_table(file, delimiter=' ', names=['Time', os.path.split(file)[1]], skiprows=1)
        df_final = df_final.merge(file_df, on='Time', how='outer')
#print(data_files)

    means = df_final.drop('Time', axis=1).mean(axis=1, skipna=True, numeric_only=True)          #calculating the mean of each row
    df_final['Mean'] = means    #adding a new column "means"
    df_final.to_csv(os.path.join(New_folder, pep_ratio + '_' + pep_name + '_' + lipid_name + '_Averages.csv'))     #saving a new file for each peptide and peptide ratio containing the normalized data and mean

    #fig = df_final.plot.scatter(x="Time", y="Mean", alpha=0.5)

    #pd.read_csv(New_folder, sep=' ', header=0)

    xvals = df_final["Time"].tolist()
    yvals = df_final["Mean"].tolist()
    #print(yvals)
    fig = plt.plot(xvals, yvals, label=filename)    #plot x,y scatter plot with lipid name

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
plt.xlabel("Time (mins)", loc='center')
plt.xlim(0, max(xmax))
plt.ylabel("Normalized Intensity", loc='center')
plt.ylim(0, 1.01)
plt.yticks(np.arange(0, 1.01, step=0.25))
plt.legend()
plt.show()
