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

def get_residue_number_and_name(df, column="P2"):
    cdata = df[column]
    resnumbers = []
    resnames = []
    for i, c in enumerate(cdata):
        print(c)
        resname = c[0]
        c2 = c[1:-1]
        resnum = int(c2)

        resnumbers.append(resnum)
        resnames.append(resname)

    df["residue_number"] = resnumbers
    df["residue_name"] = resnames
    return df


dir = "Z:\Group Share\FAM\VMD\AqpZ\MD_Simulations"

pdbfile = os.path.join(dir, "aqpz.pdb")
ppdb = pdb.PandasPdb()
ppdb.read_pdb(pdbfile)
pdbseq = ppdb.amino3to1()
pdbseq2 = ''.join(pdbseq["residue_name"])

datafile = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis\\MD_Data_MDRobin.xlsx"

lnums = [0, 1, 2, 3]
lipids = ["POCL", "TOCL", "POPE", "POPG"]
for lnum in lnums:

    pdboutput = os.path.join(dir, "AqpZ_MDRobin_" + "POCL" + str(lnum) + ".pdb")
    df = pd.read_excel(datafile)

    datacolumn = "GA"
    df = get_residue_number_and_name(df, column="P2")

    # Filter by # Lnum == 0
    df = df[df["# Lnum"] == lnum]

    resnumbers = np.array(df["residue_number"])

    for index, row in ppdb.df["ATOM"].iterrows():
        i = row["residue_number"]
        ppdb.df["ATOM"].loc[index, "b_factor"] = 0
        if i not in resnumbers:
            continue

        resname = seqdict[row["residue_name"]]
        #resname = pp.three_to_one(row["residue_name"])
        dfrow = df[df["residue_number"] == i]
        dfresname = dfrow["residue_name"]
        if np.all(resname != dfresname):
            print(resname, dfresname, i)
            continue

        try:
            ppdb.df["ATOM"].loc[index, "b_factor"] = dfrow[datacolumn].values[0]
        except:
            pass

    ppdb.to_pdb(path=pdboutput, records=None, gz=False, append_newline=True)
