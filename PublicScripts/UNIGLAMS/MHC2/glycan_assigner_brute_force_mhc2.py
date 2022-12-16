import time
import numpy as np
import matplotlib.pyplot as plt
import os
import unidec
import pandas as pd
from unidec_modules.matchtools import *


# Set the directory and engine
os.chdir("C:\\Data\\Luis Genentech\\MHC2\\")
eng = unidec.UniDec()

# Set and load the peak file
pfile = "12012022_MHC2DPA_MS2_PTCR_20221201051533_peaks2.txt"
pfile = "MHC2_DPA_Trunc_peaks.txt"
pfile = "12012022_MHC2DPA_MS2_PTCR_20221201051533_peaks.txt"
eng.load_peaks(pfile)

# Set and load the glycan file
gfile = "MHC2_DPA_Combined_tryp_glycan_list_abridged.csv"
gdf = pd.read_csv(gfile)
print(gdf.keys())

# Set the protein mass
protmass = 67987.88
protmass =  67999.04 + 3830
print("Protein Mass:", protmass)

# Set the sites
sites = ["S1", "S2", "S3", "S4"]

st = time.perf_counter()
probs_cutoff=0
# Brute force calculate all possible combinations of glycans in the list
indexes, masses, probs, glycans = get_sitematch_list(gdf, sites=sites, probs_cutoff=probs_cutoff)

# Plot Histogram
binsize = 20
bins = np.arange(np.amin(masses)+protmass, np.amax(masses)+protmass, binsize)
hist, edges = np.histogram(masses+protmass, bins=bins)
histdat = np.transpose([edges[0:-1]+binsize/2., hist])
np.savetxt(gfile[:-4]+"TheoreticalMassBruteForceHistdat"+str(probs_cutoff)+".txt", histdat)
plt.hist(masses + protmass, bins=bins)
plt.savefig(gfile[:-4]+"TheoreticalMassBruteForce"+str(probs_cutoff)+".pdf")
plt.show()


# Open the output file
outfile = gfile[:-4] + "_matches_brute_force"+str(probs_cutoff)+".xlsx"
with pd.ExcelWriter(outfile) as writer:
    # Loop through all peaks
    for mass in eng.pks.masses:
        # Calculate the difference in mass between the measured peak and the protein mass
        peakdelta = mass - protmass
        # Match the glycan masses to the delta mass
        matchindexes, matchmasses, matchprobs = sitematch_to_target(peakdelta, indexes, masses, probs, tolerance=1000)

        # Create DataFrame and load with key data
        matchdf = pd.DataFrame()
        matchdf["Measured Delta Mass"] = np.zeros_like(matchmasses) + peakdelta
        matchdf["Match Mass"] = matchmasses
        matchdf["Error"] = matchmasses - peakdelta
        matchdf["Probability"] = matchprobs
        # Normalize the probabilities
        try:
            matchdf["Sum Norm Prob"] = matchprobs / np.sum(matchprobs)
            matchdf["Max Norm Prob"] = matchprobs / np.amax(matchprobs)
        except:
            matchdf["Sum Norm Prob"] = matchprobs * 0
            matchdf["Max Norm Prob"] = matchprobs * 0
        # Find the name of the glycan for each site from the index
        for i, s in enumerate(sites):
            inds = matchindexes[:, i].astype(int)
            if len(inds) > 0:
                matchdf["Site" + str(i + 1)] = glycans[i][inds]
            else:
                matchdf["Site" + str(i + 1)] = ""

        #print(matchdf)
        # Write to excel into a sheet with the peak mass as the name
        matchdf.to_excel(writer, sheet_name=str(mass))

print("End Time:", time.perf_counter() - st)
