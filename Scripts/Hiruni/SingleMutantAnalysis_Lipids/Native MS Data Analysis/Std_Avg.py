import numpy as np
import matplotlib.pyplot as plt
import os

#Selecting the file
topdir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis2\\"

edges = ["POPG_Std", "POPE_Std", "POCL_Std"]

file = "Std_Data.xlsx"

# Calculate the average of the three standard deviation
average_std_dev = np.mean(["POPG_Std", "POPE_Std", "POCL_Std"])
