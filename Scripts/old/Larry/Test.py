import easygui
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import os
import tkinter
from tkinter import filedialog as fd
from matplotlib import rcParams

root = tkinter.Tk()
files=fd.askopenfilenames(initialdir="D:\Group Share\Fluorescence", title="Select Files", filetypes=(("txt files", "*.txt"),("csv files", "*.csv")))
files = root.tk.splitlist(files)
#print(files)

list=[]
#for f in files:
#    file = open(f)
#    lines = file.readlines()
#    file.close()

for item in files:
    file = open(item)
    lines = file.readlines()
    list=np.append(lines, skiprows=4, delimiter='\t', max_rows=len(list(files))-5)

#print(list)

for i in list:
    print(i)

