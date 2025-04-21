import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Set the directory and engine
os.chdir("C:\\Data\\Luis Genentech\\IL22_Tryp")

# Set and load the peak file
pfile = "IL22_SA4"

peaks = np.loadtxt(pfile + "_peaks.txt", delimiter=",", usecols=(0, 1))
peaks[:, 1] = peaks[:, 1] / np.sum(peaks[:, 1])

efile = "C:\Data\Luis Genentech\IL22_SA4_tryp_glycan_list_abridged_matches_brute_force_ultimate_no_cutoff.xlsx"
statsfile = pfile + "_tryp_glycan_list.csv"

sdf = pd.read_csv(statsfile)

sites = np.arange(0, 6, 2) + 1

gdata = []
pdata = []

for i, p in enumerate(peaks):

    pmass = p[0]
    pheight = p[1]
    df = pd.read_excel(efile, sheet_name=str(pmass))
    probs = df["Sum Norm Prob"].to_numpy()

    gsites = []
    psites = []
    for s in sites:
        column = "Site" + str(s)
        glycans = df[column].to_numpy()
        ug = np.unique(glycans)
        pvals = np.zeros_like(ug)
        for j, g in enumerate(ug):
            pvals[j] = np.sum(probs[glycans == g])

        gsites.append(ug)
        psites.append(pvals * pheight)

    gdata.append(gsites)
    pdata.append(psites)

unigly = np.unique(np.concatenate([np.concatenate(g) for g in gdata]))
unigly = np.concatenate([unigly, ["Unassigned"]])
print(unigly)
parray = np.zeros((len(unigly), len(sites)))
sarray = np.zeros((len(unigly), len(sites)))

# Loop through all peaks
for i, g in enumerate(gdata):
    if len(g[0]) == 0:
        for j, s in enumerate(sites):
            parray[-1, j] += peaks[i, 1]
    elif len(g) == len(sites):
        for j, s in enumerate(sites):
            pdat = pdata[i][j]
            gdat = g[j]
            for k, ug in enumerate(unigly):
                b1 = gdat == ug
                if np.any(b1):
                    parray[k, j] += np.sum(pdat[b1])
    else:
        print("ERROR:", g)

# Create a similar array from the bottom-up only
gdat = sdf["Glycan"].to_numpy()
for j, s in enumerate(sites):
    pdat = sdf["S" + str(j + 1)].to_numpy()
    for k, ug in enumerate(unigly):
        b1 = gdat == ug
        if np.any(b1):
            sarray[k, j] += np.sum(pdat[b1])

parray = parray / np.amax(parray) * 100
sarray = sarray / np.amax(sarray) * 100

plt.figure(figsize=(4.25, 5))

plt.subplot(121)
plt.imshow(parray)
ax = plt.gca()
ax.set_title("DIA-PTCR")
ax.set_xticks(range(len(sites)), labels=["1", "2", "3"])
ax.set_xlabel("Sites")
ax.set_yticks(range(len(unigly)), labels=unigly)
plt.subplot(122)
plt.imshow(sarray)
ax = plt.gca()
ax.set_title("Bottom-up")
ax.set_xticks(range(len(sites)), labels=["1", "2", "3"])
ax.set_xlabel("Sites")
ax.set_yticks(range(len(unigly)), labels=unigly)
plt.colorbar()
# plt.tight_layout()

plt.savefig(pfile + "_ultimate_barcode.pdf")
plt.savefig(pfile + "_ultimate_barcode.png")
plt.show()
