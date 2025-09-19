import pyteomics.mass as ms
import numpy as np
import os
import time
from collections import Counter
import pandas as pd
import fileinput
from isogen_tools import *
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib as mpl
import random
from Scripts.JGP.IsoGen_Analysis.unimod_obo_parser import *
from unidec.IsoDec.IsoGen.isogenc import fft_pep_seq_to_dist


atoms_to_ignore = ["Br", "I", "Cl", "F", "Hg", "Mo", "Se", "B", "Cu", "Si", "As", ]


def get_big_seq(seqs):
    big_seq = [char for seq in seqs for char in seq]
    return big_seq

def seq_to_dist_vecs_seqs(seq):
    dist = fft_pep_seq_to_dist(seq)
    vec = peptide_to_vector(seq)
    mass = peptide_to_mass(seq)

    return dist, vec, seq, mass

def parse_tsv_file(fname):
    df = pd.read_csv(fname, sep="\t")
    seqs = df["Sequence"].to_numpy()
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
    masses = []

    with multiprocessing.Pool(processes=8) as pool:
        results = pool.map(seq_to_dist_vecs_seqs, seqs)

    for result in results:
        dist, vec, seq, mass = result
        if dist is not None and vec is not None:
            dists.append(dist)
            vecs.append(vec)
            goodseqs.append(seq)
            masses.append(mass)


    dists = np.array(dists)
    vecs = np.array(vecs)

    goodseqs = np.array(goodseqs)
    if len(goodseqs) < 1:
        return None, None, None, None
    return goodseqs, dists, vecs, masses

def parse_file(fname, maxn=1000000, maxlen=200):
    if fname.endswith(".tsv"):
        seqs = parse_tsv_file(fname)

    elif fname.endswith(".mgf"):
        seqs = parse_mgf_file(fname, maxn=maxn)
    else:
        raise Exception("Unknown File Type")
    if maxn is not None:
        seqs = seqs[:maxn]

    seqs = [s for s in seqs if len(s) <= maxlen]
    print("Retreived Sequences:", len(seqs))
    goodseqs, dists, vecs, masses = seqs_to_vectors(seqs)

    #Get the filename without the extension
    filename = os.path.splitext(os.path.basename(fname))[0]
    np.savez_compressed(filename + ".npz", dists=dists, vecs=vecs, seqs=goodseqs, masses=masses)
    #np.savez_compressed("peptidedists_" + str(len(dists)) + ".npz", dists=dists, vecs=vecs, seqs=goodseqs)

def gen_random_prots(all_seqs, organism, iterations=10):
    #First concatenate all_seqs into a single list
    all_aas = []
    lengths = []
    for seq in all_seqs:
        all_aas.extend(list(seq))
        lengths.append(len(seq))

    goodseqs = np.array([])
    dists = np.array([])
    vectors = np.array([])

    for i in range(iterations):
        #Generate random protein sequences by shuffling all_aas and lengths
        # and then concatenating them
        current_all_aas = all_aas.copy()
        random_seqs = []
        random.shuffle(all_aas)
        for length in lengths:
            if length > len(all_aas):
                continue
            seq = ''.join(all_aas[:length])
            random_seqs.append(seq)
            #now remove the sequence from all_aas
            current_all_aas = current_all_aas[length:]

        curr_goodseqs, curr_dists, curr_vectors = seqs_to_vectors(random_seqs)
        if curr_goodseqs is None:
            print("No good sequences found in random sequences")
            continue

        curr_goodseqs = np.array(curr_goodseqs)
        curr_dists = np.array(curr_dists)
        curr_vectors = np.array(curr_vectors)
        if i == 0:
            goodseqs = curr_goodseqs
            dists = curr_dists
            vectors = curr_vectors
        else:
            #Add the random sequences, dists, and vectors to the goodseqs, dists, and vectors lists
            goodseqs = np.concatenate((goodseqs, curr_goodseqs))
            dists = np.concatenate((dists, curr_dists))
            vectors = np.concatenate((vectors, curr_vectors))
        print("Completed iteration", i + 1)

    #Now save the random sequences, dists, and vectors to a file
    goodseqs = np.array(goodseqs)
    dists = np.array(dists)
    vectors = np.array(vectors)
    np.savez_compressed("training_random_prots_" + str(organism) + ".npz", dists=dists, vecs=vectors, seqs=goodseqs)
    return

def mod_to_chemicalformula(mod):
    return ''.join(f"{key}{value}" for key, value in mod.composition.items())

def gen_random_seqs_even_length(all_seqs, organism, n=10000, min_length=1, max_length=200):
    big_seq = get_big_seq(all_seqs)
    np.random.shuffle(big_seq)

    random_seqs = []
    good_seqs = []
    vectors = []
    masses = []
    dists = []


    for length in range(min_length, max_length+1):
        current = 0
        index = 0

        while current < n:
            seq = ''.join(big_seq[index:index+length])

            random_seqs.append(seq)
            current += 1
            index += length

            if index + length > len(big_seq):
                np.random.shuffle(big_seq)
                index = 0

    print("Processing sequences...")
    with multiprocessing.Pool(processes=8) as pool:
        results = pool.map(seq_to_dist_vecs_seqs, random_seqs)


    print("Parsing results...")
    for r in results:
        if r[0] is not None and r[1] is not None and r[2] is not None and r[3] is not None:
            dists.append(np.array(r[0]))
            vectors.append(r[1])
            good_seqs.append(str(r[2]))
            masses.append(r[3])

    np.savez_compressed("training_random_"+ str(organism)+ "_proteins_"+str(n)+ "_min_" + str(min_length) + "_max_" +
                        str(max_length) + ".npz",
                        dists=dists, vecs=vectors, seqs=good_seqs,masses=masses)

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    os.chdir(r"C:\Users\Admin\Desktop\martylab\IsoGen\IntactProtein\Training")

    min_length = 51
    max_length = 300

    human_prot_data = np.load("human_protein_seqs.npz")
    gen_random_seqs_even_length(human_prot_data["seqs"], "human", min_length=min_length, max_length=max_length)

    yeast_prot_data = np.load("yeast_protein_seqs.npz")
    gen_random_seqs_even_length(yeast_prot_data["seqs"], "yeast", min_length=min_length, max_length=max_length)

    ecoli_prot_data = np.load("ecoli_protein_seqs.npz")
    gen_random_seqs_even_length(ecoli_prot_data["seqs"], "ecoli", min_length=min_length, max_length=max_length)


    mouse_prot_data = np.load("mouse_protein_seqs.npz")
    gen_random_seqs_even_length(mouse_prot_data["seqs"], "mouse", min_length=min_length, max_length=max_length)











