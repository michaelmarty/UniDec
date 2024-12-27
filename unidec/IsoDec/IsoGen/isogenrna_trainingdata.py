import numpy as np
import os
from isogen_tools import *

# Create all possible combinations for peptides of lengths 1 to 5
def create_rnas():
    rnaseqs = []
    for a in rnas:
        rnaseqs.append(a)

    for a in rnas:
        for b in rnas:
            rnaseqs.append(a + b)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                rnaseqs.append(a + b + c)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                for d in rnas:
                    rnaseqs.append(a + b + c + d)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                for d in rnas:
                    for e in rnas:
                        rnaseqs.append(a + b + c + d + e)

    for a in rnas:
        rnaseqs.append(a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
        )

    return rnaseqs

def create_rand_rnas(n=1000, start=6, maxlen=1000):
    seqs = []
    for i in range(n):
        length = np.random.randint(start, maxlen)
        seq = "".join(np.random.choice(["A", "C", "G", "U"], length))
        seqs.append(seq)
    return seqs

def seqs_to_vectors(seqs):
    goodseqs = []
    dists = []
    vecs = []

    for seq in seqs:
        try:
            dist = rnaseq_to_dist(seq)
            vec = rnaseq_to_vector(seq)
            goodseqs.append(seq)
            dists.append(dist)
            vecs.append(vec)
        except:
            pass

    dists = np.array(dists)
    vecs = np.array(vecs)
    goodseqs = np.array(goodseqs)
    if len(goodseqs) < 1:
        return None
    return goodseqs, dists, vecs

os.chdir("Z:\Group Share\JGP\PeptideTraining")

rnas = create_rnas()
n=100000
rnas2 = create_rand_rnas(n, maxlen=100)

# merge two lists
rnas = rnas + rnas2

print("Number of RNAs:", len(rnas))
goodseqs, dists, vecs = seqs_to_vectors(rnas)
print("Vecs Created")

np.savez_compressed(
    "rnadists_synthetic_" + str(len(dists)) + ".npz",
    dists=dists,
    vecs=vecs,
    seqs=goodseqs,
)
