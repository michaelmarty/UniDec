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
all_data = all_data[:72]
times = times[:72]

set_names = ["Blank", "0% IPTG Vector", "1% IPTG Vector", "3% IPTG Vector", "0% IPTG AqpZ_WT", "1% IPTG AqpZ_WT",
             "3% IPTG AqpZ_WT", "0% IPTG AqpZ_W14A", "1% IPTG AqpZ_W14A", "3% IPTG AqpZ_W14A", "0% IPTG AqpZ_TM",
             "1% IPTG AqpZ_TM", "3% IPTG AqpZ_TM"]

colors = ["k", "r", "orange", "salmon", "g", "lime", "y", "b", "teal", "c", "m", "orchid", "purple"]

sets = [["A1", "A2", "A3"], ["E1", "F1", "G1"], ["E2", "F2", "G2"], ["E3", "F3", "G3"], ["E4", "F4", "G4"],
        ["E5", "F5", "G5"], ["E6", "F6", "G6"], ["E7", "F7", "G7"], ["E8", "F8", "G8"],["E9", "F9", "G9"],
        ["E10", "F10", "G10"], ["E11", "F11", "G11"], ["E12", "F12", "G12"]]


set_names = ["Blank",  "1% IPTG Vector", "3% IPTG Vector", "0% IPTG AqpZ_WT", "1% IPTG AqpZ_WT",
              "0% IPTG AqpZ_W14A", "1% IPTG AqpZ_W14A", "3% IPTG AqpZ_W14A",  "3% IPTG AqpZ_TM"]

colors = ["k", "r", "orange", "salmon", "g", "lime", "y", "b", "teal", "c", "m", "orchid", "purple"]

sets = [["A1", "A2", "A3"], ["E2", "F2", "G2"], ["E3", "F3", "G3"], ["E4", "F4", "G4"],
        ["E5", "F5", "G5"], ["E7", "F7", "G7"], ["E8", "F8", "G8"],["E9", "F9", "G9"],
         ["E12", "F12", "G12"]]


fits = []

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
    print(s, j, setname)
    print(fit)
    fits.append(fit)
    plt.errorbar(xdata, ydata, yerr= ystd, label=setname, marker='o', linestyle='None', color=color)
    plt.plot(xdata, fitdata, label=setname + " Fit", color=color)
plt.legend()


plt.figure()
plt.subplot(221)
plt.title("Mid")
plt.bar(np.arange(len(fits)), [f[0] for f in fits], tick_label=set_names)

plt.subplot(222)
plt.title("Slope")
plt.bar(np.arange(len(fits)), [f[1] for f in fits], tick_label=set_names)

plt.subplot(223)
plt.title("A span")
plt.bar(np.arange(len(fits)), [f[2] for f in fits], tick_label=set_names)

plt.subplot(224)
plt.title("Background")
plt.bar(np.arange(len(fits)), [f[3] for f in fits], tick_label=set_names)
plt.tight_layout()

plt.show()
