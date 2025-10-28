import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from biopandas import pdb
from Bio import SeqIO
import Bio.PDB.Polypeptide as pp
import os, fnmatch


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


dir = "P:\Group Share\Megan E\EV_couplings\glycolysis\hexokinase_1_500_G3P_dehydrogenase"

# dir = P:\Group Share\Megan E\EV_couplings\glycolysis\hexokinase_1_500_G3P_dehydrogenase


file = find("*CouplingScores.csv", os.path.join(dir, "couplings"))[0]

pdbfile = os.path.join(dir, "1u8f_chain_O.pdb")
pdboutput = os.path.join(dir, "protein2.pdb")
fastafile = os.path.join(dir, "hexokinase_1_500_G3P_dehydrogenase.fasta")
outfile = os.path.join(dir, "probs_hexokinase_1_500_G3P_dehydrogenase.txt")

sequences = [record for record in SeqIO.parse(fastafile, "fasta")]

# Change these
s1index = 1

if s1index == 1:
    s2index = 0
    seqindex = "j"
if s1index == 0:
    seqindex = "i"
    s2index = 1

seq1 = sequences[s1index].seq
seq2 = sequences[s2index].seq
aa1 = np.arange(1, len(seq1) + 1)
aa2 = np.arange(1, len(seq2) + 1)
print(len(aa1), len(aa2))
df = pd.read_csv(file)

boo1 = df["segment_i"] != df["segment_j"]
df2 = df[boo1]

his, bins = np.histogram(df2[seqindex], bins=aa1, weights=df2["probability"])
harray = np.transpose([bins[:-1], his])

sarray = harray[harray[:, 0].argsort()]
# print(sarray)

plt.subplot(121)
plt.plot(bins[:-1], his)
plt.xlabel("Residue Number")

np.savetxt(outfile, harray)

aas1 = np.append(np.unique(np.unique(sequences[0].seq)), "Z")
his2, bins2 = np.histogram(df2["A_" + seqindex], bins=aas1, weights=df2["probability"])



plt.subplot(122)
plt.plot(aas1[:-1], his2)
plt.xlabel("Residue Name")

plt.show()

ppdb = pdb.PandasPdb()
ppdb.read_pdb(pdbfile)

pdbseq = ppdb.amino3to1()
pdbseq2 = ''.join(pdbseq["residue_name"])
starttest = 4
findlength = 10
subseq = pdbseq2[starttest:starttest + findlength]

offset = seq1.find(subseq)
if offset == -1:
    print("ERROR SEQUENCE NOT ALLIGNED")
    print(pdbseq2)
    print(seq1)
    offset = seq1.find(subseq) - starttest
    print(subseq, offset)
    print(seq1[offset:offset + findlength])
    exit()

offset = seq1.find(subseq) - starttest
print(subseq, offset)
print(seq1[offset:])
print(pdbseq2)
print(pdbseq)

for index, row in ppdb.df["ATOM"].iterrows():
    i = row["residue_number"] - 1 + offset
    try:
        resname = pp.three_to_one(row["residue_name"])
        if resname != seq1[i]:
            print(resname, seq1[i], i)
    except:
        print("Unknown Res:", resname, i, row["residue_number"])
    try:
        ppdb.df["ATOM"].loc[index, "b_factor"] = sarray[i, 1] / np.amax(sarray[:, 1])
    except:
        pass

ppdb.to_pdb(path=pdboutput, records=None, gz=False, append_newline=True)
