import requests
import os
import pandas as pd

def fetch_gpcrdb_data(code, params=None):
    base_url = "https://gpcrdb.org/services/protein/"
    url = base_url + code + "/"
    response = requests.get(url, params=params)
    response.raise_for_status()  # Raise an error for bad responses
    return response.json()

os.chdir(r"C:\Data\GPCR")

listfile = r"GPCRdb2.xlsx"
df = pd.read_excel(listfile, dtype=str)
print(df.columns)
seqs = []
accessions = []
for index, row in df.iterrows():
    # if index > 10:
    #     seqs.append("")
    #     accessions.append("")
    #     continue

    uniprot = row["GPCRs\n(UniProt)"]
    uniprot = uniprot.strip().lower() + "_human"
    print(f"Processing UniProt ID: {uniprot}")
    data = fetch_gpcrdb_data(uniprot)
    print(data)
    try:
        seq = data["sequence"]
        acc = data["accession"]
        seqs.append(seq)
        accessions.append(acc)
    except Exception as e:
        print(f"Error fetching data for {uniprot}: {e}")
        seqs.append("")
        accessions.append("")



df["Accession"] = accessions
df["Sequence"] = seqs
outfile = r"GPCRdb_with_sequences.xlsx"
df.to_excel(outfile, index=False)

df = pd.read_excel(r"GPCRdb_with_sequences.xlsx", dtype=str)
print(df.columns)

# Convert these all to Fasta
fasta_outfile = r"GPCRdb_sequences.fasta"
with open(fasta_outfile, "w") as f:
    for index, row in df.iterrows():
        uniprot = row["GPCRs\n(UniProt)"]
        acc = row["Accession"]
        seq = row["Sequence"]
        family = row["Receptor family"]
        type = row["Ligand type"]
        c = row["Class"]
        if seq == "":
            continue
        f.write(f">{acc}|{uniprot}|{family}|{type}|{c}\n")
        # Wrap sequence to 60 characters per line
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
print("Fasta file written to", fasta_outfile)
