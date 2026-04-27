import math
import multiprocessing
import os
from isogen_tools import *
from unidec.IsoDec.IsoGen.isogenc import *
from itertools import combinations_with_replacement
import random

nucleotides = ["A", "C", "G", "U"]

def get_big_seq(seqs):
    big_seq = [char for seq in seqs for char in seq]
    return big_seq

def gen_random_bigseq(a_frac, c_frac, g_frac, u_frac, len=100000):
    seq = ""
    seq += round(a_frac*len) * "A"
    seq += round(c_frac*len) * "C"
    seq += round(g_frac*len) * "G"
    seq += round(u_frac*len) * "U"

    seq = [char for char in seq]
    np.random.shuffle(seq)
    return seq

def seq_to_dist_vecs_seqs(seq):
    dist = fft_rna_seq_to_dist(seq, isolen=128)
    vec = rnaseq_to_vector(seq)
    mass = rnaseq_to_mass(seq)
    return dist, vec, seq, mass

# Create all possible combinations for RNAs of lengths 1 to 5
def create_rnas():
    print("Creating synthetic short RNAs...")
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

def create_all_rnas(min_length=1, max_length=10, max_count=500):
    seqs = []
    for i in range(min_length, max_length+1):
        print("Processing length: ", str(i))
        all_combos = []
        for combo in combinations_with_replacement(nucleotides, i):
            all_combos.append(''.join(combo))
        if len(all_combos) > max_count:
            rand_indices = random.sample(range(len(all_combos)), max_count)
            all_combos = [all_combos[i] for i in rand_indices]
        seqs.extend(all_combos)
    print("Produced", str(len(seqs)), "sequences.")
    return seqs

def create_rand_rnas(n=1000, start=6, maxlen=500):
    seqs = []
    for i in range(n):
        if i % 10000 == 0:
            print("Creating random RNA number", i)
        length = np.random.randint(start, maxlen)
        seq = "".join(np.random.choice(["A", "C", "G", "U"], length))
        seqs.append(seq)
    return seqs

def seqs_to_vectors(seqs):
    goodseqs = []
    dists = []
    vecs = []

    for i, seq in enumerate(seqs):
        if i% 10000 == 0:
            print("Processing RNA number", i)
        try:
            dist = rnaseq_to_dist(seq)
            vec = rnaseq_to_vector(seq)
            goodseqs.append(seq)
            dists.append(dist)
            vecs.append(vec)
        except:
            print("Failed to process", seq)
            pass

    dists = np.array(dists)
    vecs = np.array(vecs)
    goodseqs = np.array(goodseqs)
    if len(goodseqs) < 1:
        return None
    return goodseqs, dists, vecs

def gen_random_seqs_even_length(n=100, min_length=1, max_length=200, nuc_ratios=[0.25, 0.25, 0.25, 0.25]):
    big_seq = gen_random_bigseq(nuc_ratios[0], nuc_ratios[1], nuc_ratios[2], nuc_ratios[3])
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


    nuc_ratios_string = str(nuc_ratios[0]) + "_" + str(nuc_ratios[1]) + "_" + str(nuc_ratios[2]) + "_" + str(nuc_ratios[3])
    nuc_ratios_string = nuc_ratios_string.replace('.','-')

    np.savez_compressed("assessment_random_RNAs_"+str(n)+ "_min_" + str(min_length) + "_max_" +
                        str(max_length) + "_" + nuc_ratios_string + ".npz",
                        dists=dists, vecs=vectors, seqs=good_seqs,masses=masses)

def create_all_rnas_oflength(length):
    seqs = []
    for combo in combinations_with_replacement(nucleotides, length):
        seqs.append(''.join(combo))
    return seqs

def create_random_seq(length):
    random.shuffle(nucleotides)
    remainder = length
    seq = ''
    for i in range(len(nucleotides) - 1):
        n = random.randint(0, remainder)
        seq = seq + (n * nucleotides[i])

        remainder -= n
        if remainder == 0:
            break

    seq = seq + (remainder * nucleotides[3])
    vec = rnaseq_to_vector(seq)
    return seq, vec

def compare_vecs(v1, v2):
    #Returns 0 if they do not match and 1 if they do
    diff = v1 - v2
    for d in diff:
        if d != 0:
            return 0
    return 1


def check_exists(vec, existing_vecs):
    vec = np.array(vec)
    for v in existing_vecs:
        if compare_vecs(vec, np.array(v)):
            return True
    return False


def create_n_random_rnas(length, n=1000):
    seqs = []
    vecs = []
    #First check if there are sufficient combinations to need to not just make all of them
    combinations = math.factorial(4 + length - 1)/(math.factorial(length)*(math.factorial(4-1)))

    #in this case, generate all possible combinations
    if combinations < n:
        seqs = create_all_rnas_oflength(length)
    else:
        for i in range(n+1):
            randseq, vec = create_random_seq(length)
            if not check_exists(vec, vecs):
                seqs.append(randseq)
                vecs.append(vec)
    return seqs



if __name__ == "__main__":
    os.chdir(r"C:\Users\Admin\Documents\martylab\RNA_SeqData\Assessment")

    if False:
        nuc_ratios = [[0.25, 0.25, 0.25, 0.25],
                      [0.7, 0.1, 0.1, 0.1],
                      [0.1, 0.7, 0.1, 0.1],
                      [0.1, 0.1, 0.7, 0.1],
                      [0.1, 0.1, 0.1, 0.7]]

        n = 10
        min_length = 5
        max_length = 500

        for ratio in nuc_ratios:
            gen_random_seqs_even_length(n=n, min_length=min_length, max_length=max_length, nuc_ratios=ratio)

    if False:
        min_length = 5
        max_length = 500
        max_count = 10

        all_seqs = create_all_rnas(min_length=min_length, max_length=max_length, max_count=max_count)

        dists = []
        vectors = []
        good_seqs = []
        masses = []

        print("Processing sequences...")
        with multiprocessing.Pool(processes=8) as pool:
            results = pool.map(seq_to_dist_vecs_seqs, all_seqs)

        print("Parsing results...")
        for r in results:
            if r[0] is not None and r[1] is not None and r[2] is not None and r[3] is not None:
                dists.append(np.array(r[0]))
                vectors.append(r[1])
                good_seqs.append(str(r[2]))
                masses.append(r[3])

        np.savez_compressed("synthetic_RNAs_min" + str(min_length) + "_max" + str(max_length) + ".npz",
                            dists=dists, vecs=vectors, seqs=good_seqs, masses=masses)

    if True:
        seqs = []
        min_length = 5
        max_length = 500
        n =  500
        for i in range(min_length, max_length+1):
            print("Processing length ", str(i))
            new_seqs = create_n_random_rnas(length=i, n=n)
            seqs.extend(new_seqs)

        dists = []
        vectors = []
        good_seqs = []
        masses = []

        print("Processing sequences...")
        with multiprocessing.Pool(processes=8) as pool:
            results = pool.map(seq_to_dist_vecs_seqs, seqs)

        print("Parsing results...")
        for r in results:
            if r[0] is not None and r[1] is not None and r[2] is not None and r[3] is not None:
                dists.append(np.array(r[0]))
                vectors.append(r[1])
                good_seqs.append(str(r[2]))
                masses.append(r[3])

        np.savez_compressed("synthetic_RNAs_min" + str(min_length) + "_max" + str(max_length) + ".npz",
                            dists=dists, vecs=vectors, seqs=good_seqs, masses=masses)




