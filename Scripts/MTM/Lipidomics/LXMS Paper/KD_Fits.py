import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
import scipy.stats as stats
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

def func(x, kd, bmax):
    return bmax * x / (kd + x)


def fit_kd(x, y):
    popt, pcov = curve_fit(func, x, y, p0=[80, 500])
    return popt


path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\LX-MS\\SI_Scin Assay\\nanoSPA_Results.xlsx"

df = pd.read_excel(path)

# Remove the first row
df = df.iloc[1:]

types = ["Brain", "POPC"]
samples = ["A", "B", "C"]

x = df["Ach Conc"].to_numpy()

colors = ["k", "b", "g"]

labels = ["A", "B"]

set_names = []
set_fits = []
figure = plt.figure()
for i, t in enumerate(types):
    fits = []
    plt.subplot(1, 2, i+1)
    # Add text markers for the names
    plt.text(0.05, 0.95, labels[i], transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    for j, s in enumerate(samples):
        set_name = t + " " + s
        set_names.append(set_name)

        y = df[set_name].to_numpy()
        y = y / np.max(y)

        fit_params = fit_kd(x, y)
        fits.append(fit_params)
        #print(fit_params)
        yfit = func(x, *fit_params)

        plt.plot(x, y, marker="o", linestyle="", label=set_name, color=colors[j])
        plt.plot(x, yfit, linestyle="--", color=colors[j])
        plt.title("AChR in " + set_name[:-2] + " Nanodiscs")
        plt.xlabel("ACh Concentration (nM)")
        if i ==0:
            plt.ylabel("Normalized Fraction Bound")

    set_fits.append(fits)


set_fits = np.array(set_fits)
a = set_fits[0,:,0]
b = set_fits[1,:,0]

print(stats.ttest_ind(a, b))

means = np.mean(set_fits, axis=1)[:,0]
stds = np.std(set_fits, axis=1)[:,0]





plt.savefig(os.path.join(os.path.dirname(path), "nanoSPA_Fits.pdf"))
print(means, stds)
plt.show()