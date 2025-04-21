import time

import numpy as np
import matplotlib.pyplot as plt
import os
from unidec import engine
import pandas as pd


def random_indexes(probs, N=100):
    shape = np.shape(probs)
    indexes = np.arange(shape[1])
    rindexes = np.random.choice(indexes, p=probs, size=())


def pick_random_glycan(gdf, sites=["S1", "S1", "S2", "S2", "S3", "S3"]):
    indexes = np.arange(len(gdf))
    rindexes = []
    rmasses = []
    prob = 1
    for s in sites:
        probs = gdf[s].to_numpy() / 100.
        rindex = np.random.choice(indexes, p=probs)
        rmass = gdf["Monoisotopic mass"][rindex]
        p = probs[rindex]
        prob *= p
        rindexes.append(rindex)
        rmasses.append(rmass)
    rindexes.append(np.sum(rmasses))
    rindexes.append(prob)
    return rindexes


def sort_and_unduplicate(indexes, masses, probs, sites=["S1", "S1", "S2", "S2", "S3", "S3"]):
    # First, sort the indexes
    sindex = np.sort(indexes, axis=1)
    ### WARNING ON THIS PART. NEED To rewrite

    # Then, remove duplicates
    sunique, sunique_indexes = np.unique(sindex, return_index=True, axis=0)
    munique = masses[sunique_indexes]
    punique = probs[sunique_indexes]

    # Finally, sort everything by mass
    sortindex = np.argsort(munique)
    smasses = munique[sortindex]
    sindex = sunique[sortindex]
    sprobs = punique[sortindex]

    return sindex, smasses, sprobs


def get_random_list(gdf, sites=["S1", "S1", "S2", "S2", "S3", "S3"], N=1000):
    randoms = np.array([pick_random_glycan(gdf, sites=sites) for n in range(N)])
    masses = randoms[:, -2]
    indexes = randoms[:, :-2]
    probs = randoms[:, -1]

    #indexes, masses, probs = sort_and_unduplicate(indexes, masses, probs, sites=sites)
    return indexes, masses, probs


def match_to_random(mass, indexes, masses, probs, tolerance=5, probs_cutoff=0.01):
    b1 = np.abs(masses - mass) < tolerance
    try:
        b2 = probs / np.amax(probs) > probs_cutoff
        b1 = b1 * b2
    except: pass
    #print(np.amin(np.abs(masses - mass)), np.sum(b1), b1)
    return indexes[b1], masses[b1], probs[b1]


os.chdir("C:\\Data\\Luis Genentech")

eng = engine.UniDec()

pfile = "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_peaks.txt"
# pfile = "20220113144226FcFusionProtein_SA8_20uScans_peaks.txt"
pfile = "211029_Wendy_Isle22FC_4_warmup_PTR2_peaks.txt"

gfile = "IL22_SA4_tryp_glycan_list_abridged.csv"

gdf = pd.read_csv(gfile)
print(gdf.keys())

eng.load_peaks(pfile)

ofile = "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile.dat"
eng.load_ofile(ofile)
# print(eng.config.oligomerlist)

eng.make_oligomers(isolated=True)

protmass = np.amax(eng.olg.oligomasslist)

#deltamass = eng.pks.masses - protmass
# print(deltamass)
#peakmass = eng.pks.masses[np.argmax([p.height for p in eng.pks.peaks])]
#peakdelta = peakmass - protmass
#print(peakmass, peakdelta)

sites = ["S1", "S1", "S2", "S2", "S3", "S3"]
N = 1000000

st = time.perf_counter()
indexes, masses, probs = get_random_list(gdf, sites=sites, N=N)

plt.hist(masses+protmass, bins=N/10)
plt.savefig("TheoreticalMassMonteCarlo.pdf")
plt.show()


exit()
outfile = gfile[:-4]+"_matches.xlsx"
with pd.ExcelWriter(outfile) as writer:
    for mass in eng.pks.masses:
        peakdelta = mass - protmass
        #print(peakdelta)
        matchindexes, matchmasses, matchprobs = match_to_random(peakdelta, indexes, masses, probs)

        matchdf = pd.DataFrame()
        matchdf["Measured Delta Mass"] = np.zeros_like(matchmasses) + peakdelta
        matchdf["Match Mass"] = matchmasses
        matchdf["Error"] = matchmasses - peakdelta
        matchdf["Probability"] = matchprobs
        try:
            matchdf["Sum Norm Prob"] = matchprobs / np.sum(matchprobs)
            matchdf["Max Norm Prob"] = matchprobs / np.amax(matchprobs)
        except:
            matchdf["Sum Norm Prob"] = matchprobs * 0
            matchdf["Max Norm Prob"] = matchprobs * 0
        for i, s in enumerate(sites):
            inds = matchindexes[:, i].astype(int)
            glycans = gdf["Glycan"][inds].to_numpy()
            matchdf["Site" + str(i + 1)] = glycans

        #print(matchdf)

        matchdf.to_excel(writer, sheet_name=str(mass))

print("End Time:", time.perf_counter() - st)
