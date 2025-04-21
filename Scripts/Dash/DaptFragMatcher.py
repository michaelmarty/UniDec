import numpy as np
from DaptFragmenter import *
import unidec.tools as ud
import matplotlib.pyplot as plt
import matplotlib
from copy import deepcopy

matplotlib.use("WXAgg")


def match_to_fragments(data, fragments, tolerance=0.5, mods={"": 0}):
    matches = []
    for m in mods:
        modmass = mods[m]

        fragmasses = fragments[:, 1].astype(float) + modmass
        for d in data:
            nearindex = ud.nearest(fragmasses, d[0])
            nearmass = fragmasses[nearindex]
            error = nearmass - d[0]
            if np.abs(error) < tolerance:
                print(d[0], fragmasses[nearindex], fragments[nearindex], m)
                f = deepcopy(fragments[nearindex])
                f[0] = f[0] + m
                match = np.concatenate([f, [str(d[0]), str(d[1]), str(error)]])
                matches.append(match)
    matches = np.array(matches)
    return matches


fragments = get_all_fragments()
print(fragments)
fdata = np.transpose([fragments[:, 1].astype(float), np.ones(len(fragments))])

file = "C:\\Users\\marty\\OneDrive\\Desktop\\SYRSO230120_03-Dapto-pepsin-combine-neutral.txt"

data = np.loadtxt(file, skiprows=4)
data = data[data[:, 0] < 1700]
data[:, 1] /= np.amax(data[:, 1])
data = data[data[:, 1] > 0.10]
# print(data)

matches = match_to_fragments(data, fragments, mods=fragmods)

plt.figure(figsize=(18,9))
plt.stem(data[:, 0], data[:, 1])
plt.stem(fdata[:, 0], fdata[:, 1], markerfmt=" ", linefmt="y:")

for m in matches:
    matchmass = float(m[1])
    measuredmass = float(m[3])
    measuredint = float(m[4])
    name = m[0]
    # plt.stem(matchmass, 1, linefmt="r")
    plt.stem(measuredmass, measuredint, linefmt="r")
    plt.text(measuredmass, measuredint, name, rotation=35)

'''
plt.figure()
diffs = np.subtract.outer(data[:, 0], fdata[:, 0])
weights = np.multiply.outer(data[:, 1], fdata[:, 1])
b1 = diffs < 500
b2 = diffs > -500
b3 = b1 * b2
diffs = diffs[b3]
weights = weights[b3]
plt.hist(np.ravel(diffs), bins=1000, weights=weights)'''
manager = plt.get_current_fig_manager()
manager.frame.Maximize(True)
plt.show()
