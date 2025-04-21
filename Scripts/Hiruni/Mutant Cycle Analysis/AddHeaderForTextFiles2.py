import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import os
import pandas as pd


# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 1


# rcParams['font.size']=8

path = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_FAM_MAE_MutantCyclePaper\\SMA_4_DG_Results"
os.chdir(path)

for file in os.listdir():
    if file.endswith(".txt"):
        file_path = f"{path}\{file}"

    with open(file_path, "r+") as f:
        a = f.read()
        print(a)

        with open(file_path, "w+") as f:
            f.write("       ddG_rep1                ddG_rep2                    ddG_rep3              ddG_Average        Confidence_Interval(95%)\n" + a)
            with open(file_path, "r+") as f:
                b = f.read()
                print(b)

    #print(path)

