import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from biopandas import pdb
from Bio import SeqIO, SeqUtils
import Bio.PDB.Polypeptide as pp
import os, fnmatch

seqdict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


dir = "Z:\Group Share\FAM\VMD\AqpZ"


pdbfile = os.path.join(dir, "aqpz.pdb")
pdboutput = os.path.join(dir, "protein2.pdb")

data = [[13, "F", 0.5], [14, "W", 1], [196, "F", -.7], [224, "R", 0.25]]
df = pd.DataFrame(data, columns=["residue_number", "residue_name", "b_factor"])
print(df)


ppdb = pdb.PandasPdb()
ppdb.read_pdb(pdbfile)

pdbseq = ppdb.amino3to1()
pdbseq2 = ''.join(pdbseq["residue_name"])

#starttest=10
#subseq =pdbseq2[starttest:starttest+10]
#offset = seq1.find(subseq)-starttest
#print(subseq, offset)
#offset=0

resnumbers = np.array(df["residue_number"])
print(resnumbers)

for index, row in ppdb.df["ATOM"].iterrows():
    i = row["residue_number"]
    ppdb.df["ATOM"].loc[index, "b_factor"] = 0
    if i not in resnumbers:
        continue
    print(row["residue_name"])
    resname = seqdict[row["residue_name"]]
    #resname = pp.three_to_one(row["residue_name"])
    dfrow = df[df["residue_number"] == i]
    dfresname = dfrow["residue_name"]
    if np.all(resname != dfresname):
        print(resname, dfresname, i)
        continue

    try:
        ppdb.df["ATOM"].loc[index, "b_factor"] = dfrow["b_factor"].values[0]
    except:
        pass

ppdb.to_pdb(path=pdboutput, records=None, gz=False, append_newline=True)
