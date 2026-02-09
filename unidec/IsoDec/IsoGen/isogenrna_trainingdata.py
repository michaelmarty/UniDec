import multiprocessing
import os
from isogen_tools import *
from unidec.IsoDec.IsoGen.isogenc import *
from itertools import combinations_with_replacement

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

def create_all_rnas(min_length=1, max_length=10):
    nucleotides = ["A", "C", "G", "U"]
    seqs = []

    for i in range(min_length, max_length+1):
        for combo in combinations_with_replacement(nucleotides, i):
            seqs.append(''.join(combo))
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

def gen_random_seqs_even_length(n=100, min_length=1, max_length=200):
    big_seq = gen_random_bigseq(0.295, 0.204, 0.205, 0.296)
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


    np.savez_compressed("assessment_random_RNAs_"+str(n)+ "_min_" + str(min_length) + "_max_" +
                        str(max_length) + ".npz",
                        dists=dists, vecs=vectors, seqs=good_seqs,masses=masses)



if __name__ == "__main__":
    os.chdir(r"C:\Users\Admin\Documents\martylab\RNA_SeqData\Assessment")

    if True:
        n = 5000
        min_length = 10
        max_length = 500
        gen_random_seqs_even_length(n=n, min_length=min_length, max_length=max_length)

    if False:
        all_seqs = create_all_rnas(min_length=2, max_length=20)

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

        np.savez_compressed("synthetic_RNAs_" + str(len(good_seqs)) + ".npz",
                            dists=dists, vecs=vectors, seqs=good_seqs, masses=masses)



