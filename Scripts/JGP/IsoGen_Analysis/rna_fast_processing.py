import gzip
from Bio import SeqIO
import os
from unidec.IsoDec.IsoGen.isogen_tools import *


def add_dnaseq_to_rnadict(seq, existing_dict):
    """
    Ignore N
    Convert T to U
    Args:
        seq:
        existing_dict:

    Returns:
        dict: The updated dictionary with the new sequence added
    """
    #First replace T with U
    seq = seq.replace("T", "U")
    #Also replace N with empty string
    seq = seq.replace("N", "")
    for i, s in enumerate(seq):
        existing_dict[s] += 1

    return existing_dict

topdir = "Z:\\Group Share\\JGP\\RNA_Seq\\Data\\dna\\Chroms"

seqs = []

nucleotide_dict = {
    "A": 0,
    "C": 0,
    "G": 0,
    "U": 0,
}

#Iterate through all files in the directory
for root, dirs, files in os.walk(topdir):
    for file in files:
        print("Parsing:", file)
        if file.endswith(".gz"):
            with gzip.open(os.path.join(root, file), 'rt') as f:
                for record in SeqIO.parse(f, "fasta"):
                    nucleotide_dict = add_dnaseq_to_rnadict(str(record.seq), nucleotide_dict)

for k, v in nucleotide_dict.items():
    print(f"{k}: {v}")

total = sum(nucleotide_dict.values())

A_fract = nucleotide_dict["A"] / total
C_fract = nucleotide_dict["C"] / total
G_fract = nucleotide_dict["G"] / total
U_fract = nucleotide_dict["U"] / total

rnaveragine_comp = np.zeros(len(rna_dict_sep["A"]))

#Iterate through the keys in rna_dict_sep
for i, k in enumerate(rna_dict_sep):
    for j, v in enumerate(rna_dict_sep[k]):
        if k == "A":
            rnaveragine_comp[j] += A_fract * v
        elif k == "C":
            rnaveragine_comp[j] += C_fract * v
        elif k == "G":
            rnaveragine_comp[j] += G_fract * v
        elif k == "U":
            rnaveragine_comp[j] += U_fract * v

elements = ["C", "H", "O", "N", "P"]
print("RNA Averagine Composition:")
for i, e in enumerate(elements):
    print(f"{e}: {rnaveragine_comp[i]}")