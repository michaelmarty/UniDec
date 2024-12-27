import pyteomics.mass as ms
import numpy as np
import os
import time
from collections import Counter
import pandas as pd
import fileinput
from isogen_tools import *

os.chdir("Z:\Group Share\JGP\PeptideTraining")

def parse_tsv_file(fname):
    df = pd.read_csv(fname, sep="\t")
    seqs = df["Base Sequence"].to_numpy()
    seqs = np.unique(seqs)
    return seqs

def parse_mgf_file(fname, maxn=1000000):
    seqs = []
    count = 0
    print("Parsing:", fname)
    for line in fileinput.input(fname):
        if "SEQ=" in line:
            seq = line.split("SEQ=")[1].split("\n")[0]
            # Drop numbers and + or - values
            seq = ''.join([i for i in seq if not i.isdigit()])
            seq = seq.replace("+", "")
            seq = seq.replace("-", "")
            seq = seq.replace(".", "")

            seqs.append(seq)

            count += 1
        if count >= maxn:
            print("Max Reached:", maxn)
            break
    return np.unique(seqs)

def seqs_to_vectors(seqs):
    goodseqs = []
    dists = []
    vecs = []

    for seq in seqs:
        try:
            dist = peptide_to_dist(seq)
            vec = peptide_to_vector(seq)
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

def parse_file(fname, maxn=1000000):
    if fname.endswith(".tsv"):
        seqs = parse_tsv_file(fname)
    elif fname.endswith(".mgf"):
        seqs = parse_mgf_file(fname, maxn=maxn)
    else:
        raise Exception("Unknown File Type")
    if maxn is not None:
        seqs = seqs[:maxn]
    print("Retreived Sequences:", len(seqs))
    goodseqs, dists, vecs = seqs_to_vectors(seqs)
    np.savez_compressed("peptidedists_" + str(len(dists)) + ".npz", dists=dists, vecs=vecs, seqs=goodseqs)

if __name__ == "__main__":
    #parse_file("HeLaPeptides.tsv")
    parse_file("large_library.mgf", maxn=10000000)






