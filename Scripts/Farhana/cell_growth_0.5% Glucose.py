import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from unidec.modules import fitting

matplotlib.use('WxAgg')

file = "Z:\\Group Share\\Farhana\\Cell Growth AqpZ_16 Hr_Run2.xlsx"

# get sheet names
xl = pd.ExcelFile(file)
sheet_names = xl.sheet_names
print(sheet_names)

times = np.arange(len(sheet_names)) * 10.0

# create list of names for A-H and 1-12
names = []
for i in range(8):
    for j in range(12):
        names.append(chr(65 + i) + str(j + 1))

all_data = []

for i, s in enumerate(sheet_names):
    df = pd.read_excel(file, sheet_name=s, skiprows=3)
    # Extract row 1-8 and columns 1-12 into numpy array
    data = df.iloc[:8, 1:13].values
    fdata = np.array(data).flatten()
    all_data.append(fdata)

all_data = np.array(all_data)
all_data = all_data[:60]
times = times[:60]

set_names = ["Blank", "Vector 0.5% Glucose", "AqpZ_WT 0.5% Glucose",
             "AqpZ_W14A 0.5% Glucose", "AqpZ_TM 0.5% Glucose"]

colors = ["k", "r", "g", "b", "y"]

sets = [["A1"], ["E2", "F2", "G2"], ["E4", "F4", "G4"], ["E7", "F7", "G7"], ["E9", "F9", "G9"]]

for j, s in enumerate(sets):
    ydatas = []
    xdata = times
    setname = set_names[j]
    color = colors[j]
    for i, n in enumerate(names):
        if n not in s:
            continue
        ydata = all_data[:, i]
        ydatas.append(ydata)

    ydata = np.mean(ydatas, axis=0)
    ystd = np.std(ydatas, axis=0)
    fit, fitdata = fitting.sig_fit(xdata, ydata)
    print(fit)
    plt.errorbar(xdata, ydata, yerr= ystd, label=setname, marker='o', linestyle='None', color=color)
    plt.plot(xdata, fitdata, label=setname + " Fit", color=color)
plt.legend()
plt.show()
