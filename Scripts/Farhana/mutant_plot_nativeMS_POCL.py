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


lipid = "TOCL"

topdir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis\\" + lipid
os.chdir(topdir)



umut = ["F10A", "F13A", "F196A", "K4A", "K79A", "K155A", "R3A", "R75A", "R224A", "W14A", "W206A", "W209A", "Y100A"]
print(umut, len(umut))

# create 4x4 subplot
fig, axs = plt.subplots(4,  4, figsize=(8, 8))
for j, u in enumerate(umut):
    newdir = "WT_" + u + " Titrations All"
    os.chdir(os.path.join(topdir, newdir))
    goodfiles = [f for f in os.listdir() if f.endswith(".RAW")]
    if len(goodfiles) == 0:
        # Find directory with 25C in the name
        for d in os.listdir():
            print(newdir)
            exit()
            if "25C" in d:
                os.chdir(os.path.join(topdir, newdir, d, "rep1"))
                goodfiles = [f for f in os.listdir() if f.endswith(".RAW")]
                break
    print(newdir)
    file = goodfiles[-1]
    print(file)
    importer = UniDecImporter.get_importer(file)
    data = importer.get_avg_scan()
    importer.close()
    # find axis
    ax = axs[int(j/4), j % 4]
    ax.plot(data[:, 0], data[:, 1]/np.amax(data[:,1]), label=file, color="k")
    ax.set_xlim(6000, 10000)
    ax.set_ylim(0,1)
    ax.set_title(u)
    # ax.legend()
    if int(j/4) == 3 or j==10 or j==11 or j==9:
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
ax = axs[3,1]
ax.set_visible(False)

plt.tight_layout()

os.chdir(topdir)
plt.savefig("mutant_plot_"+lipid+".png", dpi=300)
plt.savefig("mutant_plot_"+lipid+".pdf")
plt.show()

