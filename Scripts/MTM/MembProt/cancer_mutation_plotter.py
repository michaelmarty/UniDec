import pandas as pd
import os
from pdbwriter import *

directory = "C:\\VMD\\SLCs"

os.chdir(directory)

proteins = ["AAC2", "AAC3", "UCP2"]




def get_mutation(string):
    # Example SLC25A5 R106I
    split = string.split(" ")
    mut = split[1]
    # Split mutation into residue and position
    residue1 = mut[0]
    residue2 = mut[-1]
    position = int(mut[1:-1])
    #print(residue1, position, residue2)
    return residue1, position, residue2


for p in proteins:
    file = p + "mutations.tsv"
    df = pd.read_csv(file, sep="\t")
    # Isolate type as Substitution
    sub = df[df["consequence"] == "Missense"]
    # Isolate polyphen_impact as probably_damaging
    sub = sub[sub["polyphen_impact"] == "probably_damaging"]

    outarray = []
    for s in sub["protein_change"]:
        residue1, position, residue2 = get_mutation(s)
        outarray.append([residue1, position, residue2])
    # Merge into sub
    sub["residue_name"] = [x[0] for x in outarray]
    sub["residue_number"] = [x[1] for x in outarray]
    sub["mutation_name"] = [x[2] for x in outarray]
    outdf = sub

    # Add column for histogram of residue numbers
    unique_residues = outdf["residue_number"].unique()
    residue_counts = []
    for r in unique_residues:
        count = np.sum(outdf[outdf["residue_number"] == r]["num_ssm_affected_cases"])
        residue_counts.append(count)
    residue_counts = np.array(residue_counts)
    residue_numbers = np.array(unique_residues)
    residue_names = []
    for r in residue_numbers:
        residue_names.append(outdf[outdf["residue_number"] == r]["residue_name"].values[0])
    residue_names = np.array(residue_names)
    outdf2 = pd.DataFrame({"residue_number": residue_numbers, "residue_name": residue_names, "mutation_count": residue_counts})

    print(outdf2)

    write_to_pdb(p + ".pdb", outdf2, outfile=p + "_mutations.pdb", datacolumn="mutation_count")

