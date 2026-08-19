import re
import os
import numpy as np
import matplotlib.pyplot as plt
import UniDecImporter as UniDecImporter
import matplotlib
from matplotlib import rcParams
from unidec.engine import UniDec
from UniDecImporter.ImportTools import merge_spectra

matplotlib.use('TkAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['font.size'] = 11
rcParams['font.sans-serif'] = "Arial"


lipid = "POPG"

topdir = "Z:\\Group Share\\Hiruni Jayasekera\\PAPER"
os.chdir(topdir)

umut = ["F10A", "F13A", "F196A", "K4A", "K79A", "K155A", "R3A", "R75A", "R224A", "W14A", "W206A", "W209A", "Y100A"]
print(umut, len(umut))

# create 4x4 subplot
fig, axs = plt.subplots(4,  4, figsize=(8, 8))
# Create second subplots
fig2, axs2 = plt.subplots(4,  4, figsize=(8, 8))

for j, u in enumerate(umut):
    subpath = "_" + u + "_" + lipid
    # Find files in the directory that have the subpath in the name
    goodfiles = [f for f in os.listdir() if subpath in f and f.endswith(".RAW")]
    print(u, goodfiles)
    if len(goodfiles) == 0:
        print("No files found for", u)
        raise FileNotFoundError
    elif len(goodfiles) == 3:
        file = goodfiles[0]
        importer = UniDecImporter.get_importer(file)
        data = importer.get_avg_scan()
    elif len(goodfiles) > 10:
        print("Combining Files")
        # Find all files with Rep1 in name
        rep1files = [f for f in goodfiles if "Rep2" in f]
        data_import = []
        for f in rep1files:
            importer = UniDecImporter.get_importer(f)
            data = importer.get_avg_scan()
            data_import.append(data)
            importer.close()
        data = merge_spectra(data_import)

    print(len(data))

    eng=UniDec()
    eng.pass_data_in(data, silent=False)
    eng.process_data()
    eng.run_unidec()

    decondat = eng.data.massdat

    # find axis
    ax = axs2[int(j/4), j % 4]
    ax.plot(decondat[:, 0], decondat[:, 1]/np.amax(decondat[:,1]), color="k")
    ax.set_xlim(97000, 107000)
    ax.set_ylim(0,1)
    ax.set_title(u)
    # ax.legend()
    if int(j/4) == 3 or j==10 or j==11 or j==9:
        # Set this to italic
        ax.set_xlabel("Mass (Da)")
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

    # find axis
    ax = axs[int(j / 4), j % 4]
    ax.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]), color="k")
    ax.set_xlim(6000, 10000)
    ax.set_ylim(0, 1)
    ax.set_title(u)
    # ax.legend()
    if int(j / 4) == 3 or j == 10 or j == 11 or j == 9:
        # Set this to italic
        ax.set_xlabel(r"$m/z$")
    else:
        ax.set_xticklabels([])
    ax.set_yticks([0, 0.5, 1])
    if j % 4 == 0:
        # ax.set_ylabel("Intensity")
        ax.set_yticklabels(["0", "%", "100"])
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

ax = axs2[3,2]
ax.set_visible(False)
ax = axs2[3,3]
ax.set_visible(False)
ax = axs2[3,1]
ax.set_visible(False)

fig.tight_layout()
fig2.tight_layout()

os.chdir(topdir)
fig2.savefig("mutant_plot_"+lipid+"_decon.png", dpi=300)
fig2.savefig("mutant_plot_"+lipid+"_decon.pdf")


fig.savefig("mutant_plot_"+lipid+".png", dpi=300)
fig.savefig("mutant_plot_"+lipid+".pdf")

plt.show()

