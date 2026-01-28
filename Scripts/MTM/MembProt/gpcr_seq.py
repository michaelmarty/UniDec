import os
from Bio import SeqIO
import gget
import pandas as pd

def get_from_tdf(tdf, input_col, output_col, input_value, input_mixed=False):
    if input_mixed:
        # Find rows where input_col contains input_value as substring
        result = tdf.loc[tdf[input_col].str.contains(input_value, na=False), output_col].values
    else:
        result = tdf.loc[tdf[input_col] == input_value, output_col].values
    if len(result) == 1:
        return result[0]
    elif len(result) > 1:
        print(f"Multiple entries found for {input_value} in column {input_col}. Returning first match.")
        return result[0]
    else:
        print(f"No entry found for {input_value} in column {input_col}.")
        return None

os.chdir(r"C:\Data\GPCR")

fasta_file = r"gpcr-gross-alignment.fasta"

translator_file=r"HUMAN_9606_idmapping_selected.tsv"
tdf = pd.read_table(translator_file, delimiter="\t", dtype=str)
print(tdf.columns)



for record in SeqIO.parse(fasta_file, "fasta"):
    print(record.id)
    uniprot = record.id.split("|")[1]
    print("Uniprot ID:", uniprot)
    if "." in uniprot:
        pdbid = uniprot.split(".")[0]
        print("PDB ID:", pdbid)
        # Find the Ensembl for this PDB ID
        ensembl = get_from_tdf(tdf, "PDB", "Ensembl", pdbid, input_mixed=True)
        print("Ensembl ID:", ensembl)
    else:
        ensembl = get_from_tdf(tdf, "UniProtKB-AC", "Ensembl", uniprot)

    if ensembl is None:
        continue
    else:
        # Strip after .
        if "." in ensembl:
            ensembl = ensembl.split(".")[0]

    info = gget.info(ensembl)
    print(info.keys())
    for k in info.keys():
        print(f"{k}: {info[k]}")
    exit()