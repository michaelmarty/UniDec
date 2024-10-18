import numpy as np
import os
import time
from isogen_tools import get_dist_from_formula

def parse_formulas(formulas, isolen=128, cutoff=0.001):
    dists = []
    goodforms = []
    for f in formulas:
        try:
            dist = get_dist_from_formula(f, isolen, cutoff)
            dists.append(dist)
            goodforms.append(f)
        except:
            pass

    dists = np.array(dists)
    goodforms = np.array(goodforms)
    return dists, goodforms

def cleanup_formulas(formulas):
    formulas = np.unique(formulas)
    print(len(formulas))
    # Rmove everything after + or -
    for i in range(len(formulas)):
        f = formulas[i]
        if "+" in f:
            formulas[i] = f.split("+")[0]
        if "-" in f:
            formulas[i] = f.split("-")[0]

    formulas = np.unique(formulas)
    return formulas

if __name__ == "__main__":

    starttime= time.perf_counter()

    # Search pubchem and get all atomic formulas
    os.chdir("Z:\Group Share\JGP\PubChem")
    fname = "CID-Mass.txt"
    formulas = np.genfromtxt(fname, dtype=str, delimiter="\t", max_rows=100000000, usecols=[1])
    formulas = cleanup_formulas(formulas)

    print(len(formulas))

    dists, goodforms = parse_formulas(formulas)

    np.savez_compressed("isodists_"+str(len(dists))+".npz", dists=dists, formulas=goodforms)

    print("Time:", time.perf_counter() - starttime)

