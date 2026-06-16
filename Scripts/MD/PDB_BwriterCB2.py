import pandas as pd
import numpy as np
from biopandas import pdb
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


dir = r"C:\VMD\CB2"

pdbfile = os.path.join(dir, "cb2_fullaf_aligned_1_312.pdb")
ppdb = pdb.PandasPdb()
ppdb.read_pdb(pdbfile)
pdbseq = ppdb.amino3to1()
pdbseq2 = ''.join(pdbseq["residue_name"])

lname = "POPG"

datafile = f"C:\\VMD\\CB2\\CG\\POPG5\\gromacs\\Interaction_{lname}\\Dataset_{lname}\\Dataset.csv"
# datafile = f"C:\\VMD\\CB2\\AA\\POPG5\\namd\\Interaction_{lname}\\Dataset_{lname}\\Dataset.csv"


pdboutput = os.path.join(dir, f"CB2_{lname}_cg_res.pdb")
df = pd.read_csv(datafile)

datacolumn = "Residence Time"
# datacolumn = "Occupancy"
df = get_residue_number_and_name(df, column="Residue", mode=1)

# Normalize datacolumn
array = df[datacolumn]

if datacolumn == "Residence Time":
    array /= np.amax(array)
else:
    array /= 100
df[datacolumn] = array

resnumbers = np.array(df["residue_number"])

chains = ["P"]

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
    # Get the chain id for the atom in the pdb
    pdbchain = row["chain_id"]

    # Get the chain id for the data by assuming that the chain id is in order
    for i, row in enumerate(dfrow.iterrows()):
        datachain = chains[i]
        if pdbchain == datachain:
            try:
                ppdb.df["ATOM"].loc[index, "b_factor"] = dfrow[datacolumn].values[i]
            except:
                pass
print("Writing to :", pdboutput)
ppdb.to_pdb(path=pdboutput, records=None, gz=False, append_newline=True)
