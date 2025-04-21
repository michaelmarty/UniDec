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
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "HSD": "H"}

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def get_residue_number_and_name(df, column="Residue", mode=0):
    cdata = df[column]
    resnumbers = []
    resnames = []
    for i, c in enumerate(cdata):
        print(c)
        if mode == 0:
            resname = c[0]
            c2 = c[1:-1]
            resnum = int(c2)
        if mode == 1:
            resname = c[-3:]
            resname = seqdict[resname]
            c2 = c[:-3]
            resnum = int(c2)
        else:
            print("Not recognized")
            exit()


        resnumbers.append(resnum)
        resnames.append(resname)

    df["residue_number"] = resnumbers
    df["residue_name"] = resnames
    return df


dir = "Z:\Group Share\Farhana\AAC2\MD"

pdbfile = os.path.join(dir, "ucp2fixed.pdb")
ppdb = pdb.PandasPdb()
ppdb.read_pdb(pdbfile)
pdbseq = ppdb.amino3to1()
pdbseq2 = ''.join(pdbseq["residue_name"])

datafile = "Z:\\Group Share\\MDJ\\UCP2\\MD\\UCP2_PEPCCL_CG_1\\gromacs\\Interaction_POPE\\Dataset.xlsx"


pdboutput = os.path.join(dir, "POPE_Occ" + ".pdb")
df = pd.read_excel(datafile)

datacolumn = "Occupancy"
df = get_residue_number_and_name(df, column="Residue", mode=1)

# Normalize datacolumn
array = df[datacolumn]
array /= 100
df[datacolumn] = array


resnumbers = np.array(df["residue_number"])

chains = ["A", "B", "C", "D"]

# Loop through every atom in the pdb
for index, row in ppdb.df["ATOM"].iterrows():
    # Set default b-factor to 0
    ppdb.df["ATOM"].loc[index, "b_factor"] = 0
    # Get the pdb residue number
    i = row["residue_number"]
    # If the residue number is not in the data, skip
    if i not in resnumbers:
        continue

    # Get the residue name from the sequence dictionary for the atom in the pdb
    resname = seqdict[row["residue_name"]]
    #resname = pp.three_to_one(row["residue_name"])

    # Get the row in the data that corresponds to the residue number
    dfrow = df[df["residue_number"] == i]
    dfresname = dfrow["residue_name"]
    # If the residue name in the pdb does not match the residue name in the data, print the residue name and continue
    if np.all(resname != dfresname):
        print("Sequence Resname:", resname)
        print("Res Number:", i)
        print(dfrow)
        exit()
        continue


    try:
        ppdb.df["ATOM"].loc[index, "b_factor"] = dfrow[datacolumn].values[0]
    except:
        pass

ppdb.to_pdb(path=pdboutput, records=None, gz=False, append_newline=True)

