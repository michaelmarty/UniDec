import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import os

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 1


# rcParams['font.size']=8

path = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_FAM_MAE_MutantCyclePaper\\DMC_3_Extract_grid_Extract6 Data - Copy"
os.chdir(path)

def read_text_file(file_path):
    with open(file_path, 'r') as f:


        for file in os.listdir():
            if file.endswith(".txt"):
                file_path = f"{path}\{file}"
                read_text_file(file_path)
                print(file_path)
                lines = f.readlines()
            with open(file_path, 'w') as f:
                for line in lines:
                    if line.strip("\n") != "header ":
                        f.write(line)

