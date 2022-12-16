import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Set the directory and engine
os.chdir("D:\\Data\\Luis Genentech\\IL22_Tryp")

# Set and load the peak file
pfile = "IL22_SA8"

peaks = np.loadtxt(pfile + "_peaks.txt", delimiter=",", usecols=(0, 1))
peaks[:, 1] = peaks[:, 1] / np.sum(peaks[:, 1])

efile = pfile + "_tryp_glycan_list_4thsite_matches_brute_force3.xlsx"
statsfile = pfile + "_tryp_glycan_list_4thsite.csv"

sdf = pd.read_csv(statsfile)

sites = np.arange(0, 8, 2) + 1

gdata = []
pdata = []
equaldata = []

for i, p in enumerate(peaks):

    pmass = p[0]
    pheight = p[1]
    df = pd.read_excel(efile, sheet_name=str(pmass))
    probs = df["Sum Norm Prob"].to_numpy()
    equalprobs = np.ones_like(probs) / len(probs)

    gsites = []
    psites = []
    psites2 = []
    for s in sites:
        column = "Site" + str(s)
        glycans = df[column].to_numpy()
        ug = np.unique(glycans)
        pvals = np.zeros_like(ug)
        pvals2 = np.zeros_like(ug)
        for j, g in enumerate(ug):
            pvals[j] = np.sum(probs[glycans == g])
            pvals2[j] = np.sum(equalprobs[glycans == g])

        gsites.append(ug)
        psites.append(pvals * pheight)
        psites2.append(pvals2 * pheight)

    gdata.append(gsites)
    pdata.append(psites)
    equaldata.append(psites2)

unigly = np.unique(np.concatenate([np.concatenate(g) for g in gdata]))
unigly = np.concatenate([unigly, ["Unassigned"]])
print(unigly)
parray = np.zeros((len(unigly), len(sites)))
sarray = np.zeros((len(unigly), len(sites)))
earray = np.zeros((len(unigly), len(sites)))

# Loop through all peaks
for i, g in enumerate(gdata):
    if len(g[0]) == 0:
        for j, s in enumerate(sites):
            parray[-1, j] += peaks[i, 1]
            earray[-1, j] += peaks[i, 1]
    elif len(g) == len(sites):
        for j, s in enumerate(sites):
            pdat = pdata[i][j]
            pdat2 = equaldata[i][j]
            gdat = g[j]
            for k, ug in enumerate(unigly):
                b1 = gdat == ug
                if np.any(b1):
                    parray[k, j] += np.sum(pdat[b1])
                    earray[k, j] += np.sum(pdat2[b1])
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
earray = earray / np.amax(earray) * 100

plt.figure(figsize=(7, 5))

plt.subplot(131)
plt.imshow(parray)
ax = plt.gca()
ax.set_title("UNIGLAMS")
ax.set_xticks(range(len(sites)), labels=["1", "2", "3", "4"])
ax.set_xlabel("Sites")
ax.set_yticks(range(len(unigly)), labels=unigly)
plt.subplot(132)
plt.imshow(sarray)
ax = plt.gca()
ax.set_title("Bottom-up")
ax.set_xticks(range(len(sites)), labels=["1", "2", "3", "4"])
ax.set_xlabel("Sites")
ax.set_yticks(range(len(unigly)), labels=unigly)
#plt.colorbar()
plt.subplot(133)
plt.imshow(earray)
ax = plt.gca()
ax.set_title("UNIGLAMS-Only")
ax.set_xticks(range(len(sites)), labels=["1", "2", "3", "4"])
ax.set_xlabel("Sites")
ax.set_yticks(range(len(unigly)), labels=unigly)
plt.colorbar()
# plt.tight_layout()

df2 = pd.DataFrame(sarray.transpose(), columns=unigly)
df2.to_excel(pfile+"_stats_array_barcode3.xlsx")

df2 = pd.DataFrame(parray.transpose(), columns=unigly)
df2.to_excel(pfile+"_prob_array_barcode3.xlsx")

df2 = pd.DataFrame(earray.transpose(), columns=unigly)
df2.to_excel(pfile+"_equal_array_barcode3.xlsx")

#np.savetxt(pfile+"_stats_array_barcode3.txt", sarray)
#np.savetxt(pfile+"_prob_array_barcode3.txt", parray)

plt.savefig(pfile + "_barcode_triple_3.pdf")
plt.savefig(pfile + "_barcode_triple_3.png")
plt.show()
