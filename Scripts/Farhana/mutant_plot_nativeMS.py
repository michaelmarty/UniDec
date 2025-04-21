import re
import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.UniDecImporter as UniDecImporter
import matplotlib
from matplotlib import rcParams

matplotlib.use('TkAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['font.size'] = 11
rcParams['font.sans-serif'] = "Arial"

def get_mutant_name(path):
    """
    Get the mutant name from the file path.
    :param path: File path
    :return: Mutant name (str)
    """
    name = os.path.basename(path)
    try:
        # Regex to search for pattern _ followed by a capital letter and a number and then an A and then either _ or .
        mutant = re.search(r"_([A-Z]\d+A[_.])", name).group()
        # Drop the first and last character
        mutant = mutant[1:-1]
    except:
        if "WT" in name:
            mutant = "WT"
        else:
            mutant = "Unknown"
            raise ValueError("Could not find mutant name in file path.", f)
    return mutant

topdir = "Z:\\Group Share\\Farhana\\AqpZ Base\\AqpZ\\Best"
os.chdir(topdir)

files = os.listdir(topdir)
goodfiles = []

mutants = []

for f in files:
    if not f.endswith(".RAW"):
        continue
    mutant = get_mutant_name(f)
    mutants.append(mutant)
    goodfiles.append(f)

umut = np.unique(mutants)
# Put WT first
umut = np.concatenate([["WT"], umut[umut != "WT"]])
print(umut, len(umut))
# create 4x4 subplot
fig, axs = plt.subplots(4,  4, figsize=(8, 8))
for j, u in enumerate(umut):
    for i, m in enumerate(mutants):
        if u == m:
            file = goodfiles[i]
            print(file)
            importer = UniDecImporter.get_importer(file)
            data = importer.get_avg_scan()
            importer.close()
            # find axis
            ax = axs[int(j/4), j % 4]
            ax.plot(data[:, 0], data[:, 1]/np.amax(data[:,1]), label=file, color="k")
            ax.set_xlim(5000, 10000)
            ax.set_ylim(0,1)
    ax.set_title(u)
    # ax.legend()
    if int(j/4) == 3 or j==10 or j==11:
        # Set this to italic
        ax.set_xlabel(r"$m/z$")
    else:
        ax.set_xticklabels([])
    ax.set_yticks([0, 0.5, 1])
    if j % 4 == 0:
        # ax.set_ylabel("Intensity")
        ax.set_yticklabels(["0", "%","100"])
    else:
        ax.set_yticklabels([])

    # Turn off top and right axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

ax = axs[3,2]
ax.set_visible(False)
ax = axs[3,3]
ax.set_visible(False)

plt.tight_layout()
plt.savefig("mutant_plot.png", dpi=300)
plt.savefig("mutant_plot.pdf")
plt.show()

