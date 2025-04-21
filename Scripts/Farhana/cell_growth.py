import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from unidec.modules import fitting

matplotlib.use('WxAgg')

file = "Z:\\Group Share\\Farhana\\Cell Growth AqpZ_12 Hr.xlsx"

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

set_names = ["Blank", "0% Vector", "1% Vector", "3% Vector", "0% AqpZ_WT", "1% AqpZ_WT", "3% AqpZ_WT", "0% AqpZ_W14A",
             "1% AqpZ_W14A", "3% AqpZ_W14A", "0% AqpZ_TM", "1% AqpZ_TM", "3% AqpZ_TM"]

colors = ["k", "r", "firebrick", "red", "g", "lime", "y", "b", "teal", "c", "m", "purple", "orchid"]

sets = [["A1", "A2", "A3"], ["B1", "C1", "D1"], ["B2", "C2", "D2"], ["B3", "C3", "D3"],
        ["B4", "C4", "D4"], ["B5", "C5", "D5"], ["B6", "C6", "D6"], ["B7", "C7", "D7"],
        ["B8", "C8", "D8"], ["B9", "C9", "D9"], ["B10", "C10", "D10"], ["B11", "C11", "D11"], ["B12", "C12", "D12"]]

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
