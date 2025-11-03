import numpy as np
import os
import pandas as pd

def df_to_df(df):
    # Rename Alignment ID to ID
    df = df.rename(columns={"Alignment ID": "ID"})
    # Rename "Metabolite name" to "Name"
    df = df.rename(columns={"Metabolite name": "Name"})
    # Rename Ontology to CompoundClass
    df = df.rename(columns={"Ontology": "CompoundClass"})
    # Rename Average Rt(min) to RT
    df = df.rename(columns={"Average Rt(min)": "RT"})
    # Rename "INCHIKEY" to InChIKey
    df = df.rename(columns={"INCHIKEY": "InChIKey"})
    # Keep only ID, Name, CompoundClass, RT, Reference RT, InChIKey and SMILES
    df = df[["ID", "Name", "CompoundClass", "InChIKey", "SMILES", "RT", "Reference RT"]]

    return df

def msp_to_df(mspfile):
    records = []
    with open(mspfile, 'r') as f:
        lines = f.readlines()

    record = {}
    for line in lines:
        line = line.strip()
        if line == "":
            if record:
                records.append(record)
                record = {}
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            record[key.strip()] = value.strip()
    if record:
        records.append(record)

    df = pd.DataFrame(records)
    # Keep only relevant columns and rename
    df = df.rename(columns={"NAME": "Name", "INCHIKEY": "InChIKey", "SMILES": "SMILES", "RETENTIONTIME": "RT", "COMPOUNDCLASS": "CompoundClass"})
    df = df[["Name", "CompoundClass", "InChIKey", "SMILES", "RT"]]
    return df

def get_rt_for_name(df, name, inchi, smiles, newcol="RT"):
    rt = df.loc[df['Name'] == name, newcol].values
    if len(rt) == 1:
        return rt[0]
    elif len(rt) > 1:
        # print(f"Duplicate names found for {name}, using first occurrence.")
        # Find the one with the same InChIKey or SMILES if possible
        filtered = df[(df['Name'] == name) & ((df['InChIKey'] == inchi) | (df['SMILES'] == smiles))]
        if len(filtered) == 1:
            return filtered.iloc[0][newcol]
        else:
            print(f"Could not resolve duplicates for {name} using InChIKey or SMILES, using first occurrence.")
            return rt[0]  # or handle duplicates as needed
    else:
        return None

def df_to_msp(df, mspfile, outfile, newcol="RT"):
    with open(mspfile, 'r') as f:
        lines = f.readlines()

    with open(outfile, 'w') as f:
        name = ""
        inchi = ""
        smiles = ""
        for line in lines:
            if name == "":
                f.write(line)
                if "NAME:" in line:
                    name = line.split(":", 1)[1].strip()
            else:
                if "INCHIKEY:" in line:
                    inchi = line.split(":", 1)[1].strip()
                    f.write(line)

                if "SMILES:" in line:
                    smiles = line.split(":", 1)[1].strip()
                    f.write(line)

                if "RETENTIONTIME:" in line:
                    newrt = get_rt_for_name(df, name, inchi, smiles, newcol)
                    try:
                        f.write(f"RETENTIONTIME: {float(newrt)}\n")
                    except:
                        print(f"Could not convert RT to float for {name}: {newrt}")
                        f.write(line)
                    name = ""
                    inchi = ""
                    smiles = ""
                else:
                    f.write(line)


os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")
# posdf = msp_to_df("MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp")
# # Write to CSV
# posdf.to_csv("Pos_C1_Ref_FromMSP.csv", index=False)
#
# # Read neg file
# negdf = msp_to_df("MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp")
# # Write to CSV
# negdf.to_csv("Neg_C1_Ref_FromMSP.csv", index=False)

# Test the df to msp function
posdf = pd.read_csv("Pos_C1_Ref_FromMSP.csv")
df_to_msp(posdf, "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp", "Test_Pos_C1.msp")

exit()

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

posfile = "PosIDsC1_Ref.csv"
negfile = "NegIDsC1_Ref.csv"

posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)
# Combine dataframes
df = pd.concat([negdf, posdf])

df = df_to_df(df)

outfile = "C1_Ref_Training.csv"
df.to_csv(outfile, index=False)

