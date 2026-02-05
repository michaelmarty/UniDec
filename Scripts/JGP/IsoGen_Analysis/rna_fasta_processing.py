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
    bad_nucs_dict = {}

    #First replace T with U
    seq = seq.replace("T", "U")

    for i, s in enumerate(seq):
        if s in existing_dict:
            existing_dict[s] += 1
        else:
            if s in bad_nucs_dict:
                bad_nucs_dict[s] += 1
            else:
                bad_nucs_dict[s] = 1

    for k,v in bad_nucs_dict.items():
        print("nuc: ", k, " frac: ", v/len(seq))

    return existing_dict

topdir = r"C:\Users\Admin\Documents\martylab\RNA_SeqData\HomoSapiens\Chrom_fastas"
os.chdir(topdir)

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
                    #seqs.append(str(record.seq.replace("T", "U")))
                    nucleotide_dict = add_dnaseq_to_rnadict(str(record.seq), nucleotide_dict)


for k, v in nucleotide_dict.items():
    print(f"{k}: {v}")

total = sum(nucleotide_dict.values())

A_fract = nucleotide_dict["A"] / total
C_fract = nucleotide_dict["C"] / total
G_fract = nucleotide_dict["G"] / total
U_fract = nucleotide_dict["U"] / total

print("A Frac:", A_fract)
print("C Frac:", C_fract)
print("G Frac:", G_fract)
print("U Frac:", U_fract)

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

rnaveragine_comp = np.append(rnaveragine_comp, [1])
rnaveragine_comp[3] += 4
elements = ["C", "H", "O", "N", "P"]
print("RNA Averagine Composition:")
for i, e in enumerate(elements):
    print(f"{e}: {rnaveragine_comp[i]}")