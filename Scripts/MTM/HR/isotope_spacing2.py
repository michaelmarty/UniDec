import molmass
import numpy as np
from unidec.modules.isotopetools import makemass
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.use("wxAgg")

file = "data.txt"
formulas = []
datas = []
data = []

def formula_to_avg_spacing(formula):
    formula = molmass.Formula(formula)
    spectrum = formula.spectrum().dataframe()
    diffs = []
    probs = []

    masses = spectrum["Relative mass"].tolist()
    intensities = spectrum["Intensity %"].tolist()
    for i in range(len(spectrum) - 1):
        diff = masses[i + 1] - masses[i]
        prob = intensities[i]
        diffs.append(diff)
        probs.append(prob)

    weighted_avg = sum(diffs[i] * probs[i] for i in range(len(diffs))) / sum(probs)
    return weighted_avg

def avg_spacing_from_data(data):
    data = np.array(data)
    diffs = []
    probs = []

    newdata = []

    startint = int(np.round(np.amin(data[:,0])))
    endint = int(np.round(np.amax(data[:,0])))
    startingdecimal = data[0, 0] - startint

    for i in range(int(startint), int(endint) + 1):
        b1 = data[:, 0] < startingdecimal+i + 0.5
        b2 = data[:, 0] > startingdecimal+i - 0.5
        b1 = np.logical_and(b1, b2)
        intensity = 0
        weighted_average = i
        if np.any(b1):
            intensity = np.sum(data[b1, 1])
        if intensity != 0:
            weighted_average = np.sum(data[b1, 0] * data[b1, 1]) / intensity
        newdata.append([weighted_average, intensity])
    newdata= np.array(newdata)
    data = newdata

    for i in range(len(data) - 1):
        diff = data[i + 1][0] - data[i][0]
        prob = data[i][1]
        diffs.append(diff)
        probs.append(prob)

    weighted_avg = sum(diffs[i] * probs[i] for i in range(len(diffs))) / sum(probs)

    return weighted_avg

with open(file, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

            # if line starts with "C" grab it as a formula
        if line.startswith("C"):
            formulas.append(line)
            if len(data) > 0:
                datas.append(data)
            data = []

        lines = line.split()
        if len(lines) < 2:
            continue
        try:
            lines = [float(x) for x in lines]
        except ValueError:
            continue
        if len(lines) != 2:
            continue

        data.append(lines)
datas.append(data)
print(formulas)

avg_spacings = [formula_to_avg_spacing(f) for f in formulas]
print("Average spacings for formulas:", np.array(avg_spacings))

avg_spacings2 = [avg_spacing_from_data(d) for d in datas]
print("Average spacings from data:", np.array(avg_spacings2))

x = np.arange(2, 5, 0.1)
masses = 10**x
formulas = [makemass(mass)[0] for mass in masses]

avg_spacings3 = [formula_to_avg_spacing(f)-1 for f in formulas]

plt.plot(masses, avg_spacings3, label="From formula")
plt.show()


