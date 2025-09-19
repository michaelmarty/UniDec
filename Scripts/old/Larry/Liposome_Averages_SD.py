import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter import filedialog as fd
import pandas as pd
from matplotlib import rcParams
import seaborn as sns

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
files = root.tk.splitlist(files)
root.destroy()

newpath=()
xmax=()
filename=()
pep_name=()
pep_ratio=()
lipid_name=()
data_sets=[] # Variable to save all the data sets

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
    newfile = 'Normalized' + '_' + filename[2] + '_' + filename[3] + '_' + filename[1] + '_' + filename[0] + '.txt'   #creating a new file using parts of original file name
    filepath = os.path.join(New_folder, newfile).replace("/","\\")  #Forcing all of the slashes in directory to be same direction
    newfile = open(filepath, "w+")  #opening the new file for writing
    np.savetxt(filepath, data, header='Time Intensity', comments='')     #saving normalized ydata to new file
    data_sets.append(pep_ratio + '_' + pep_name + '_' + lipid_name)

data_sets = list(np.unique(np.array(data_sets))) # Get the unique list of data sets (pep_ratio + pep_name)

df_final = pd.DataFrame(columns=['Time'])
ups= []
downs = []
for dset in data_sets:
    data_files = glob.glob(os.path.join(New_folder, 'Normalized_' + dset + '*.txt'))  # list of all the normalized data files for a given data_set
    df_set = pd.DataFrame(columns=['Time'])  # creating a pandas data frame to store the results
    print(data_files)
    for file in data_files:
        file_df = pd.read_table(file, delimiter=' ', names=['Time', os.path.split(file)[1]], skiprows=1)
        df_set = df_set.merge(file_df, on='Time', how='outer')


    means = df_set.drop('Time', axis=1).mean(axis=1, skipna=True)
    stds = df_set.drop('Time', axis=1).std(axis=1, skipna=True)
    name_col = 'Mean_' + dset
    df_set[name_col] = means
    name_col = 'Std_' + dset
    df_set[name_col] = stds
    #name_col = 'UP_' + dset
    #df_set[name_col] = means+stds
    #name_col = 'DN_' + dset
    #df_set[name_col] = means - stds
    #ups.append(means+stds)
    #downs.append(means+stds)
    df_set.to_csv(os.path.join(New_folder, dset + '_Averages.csv'), index=False)
    df_set = df_set[['Time', name_col]]  # Subsetting only columns needed for plotting
    df_final = df_final.merge(df_set, on='Time', how='outer')

df_final.to_csv('All_data_sets_averages.csv')

# transforming data from wide to long format for plotting
df_melt = df_final.melt(id_vars='Time', var_name='data_set', value_name='mean')

#df = pd.read_csv('All_data_sets_averages.csv')
#array = pd.DataFrame(pd.read_csv('All_data_sets_averages.csv', usecols=[2,3,4])).to_numpy()
#print(array)

# using sns (seaborn) because it is easier to color by different datasets
fig = sns.lineplot(data=df_melt, x='Time', y='mean', hue='data_set', alpha=0.5)
#fig = sns.lineplot(data=df_melt, x='Time', y='UP_', hue='data_set', alpha=0.5)
#fig = plt.plot(data=df_melt, x='Time', y='mean', label=filename)    #plot x,y scatter plot with lipid name

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
sns.despine(left=False, right=True, bottom=False, top=True)
plt.xlabel("Time (mins)", loc='center')
plt.xlim(0, max(xmax))
plt.ylabel("Normalized Intensity", loc='center')
plt.ylim(0, 1.01)
plt.yticks(np.arange(0, 1.01, step=0.25))
##fig.legend#(loc='center left', bbox_to_anchor=(1, 0.5)) # putting legend outside of plot
plt.show()
