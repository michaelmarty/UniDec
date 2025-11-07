import numpy as np
import os
import pandas as pd
import Scripts.MTM.Lipidomics.LipidFunctions as LipidFunctions

def df_to_df(df):
    # Clean up duplicate names with same name and INCHIKEY
    # Keep the one with highest S/N average
    print(len(df))
    df = df.sort_values(by=['S/N average'], ascending=False)
    df = df.drop_duplicates(subset=['Metabolite name', 'INCHIKEY'], keep='first')
    print(len(df))

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
    df = df[["ID", "Name", "CompoundClass", "InChIKey", "SMILES", "RT", "Reference RT", "S/N average"]] #"Tail Length", "Tail Unsaturation",

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

def get_rt_for_name(df, name, inchi, smiles, newcol="RT", verbose=True):
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
            if verbose:
                print(f"Couldn't resolve duplicates for {name}. Using first")
            return rt[0]  # or handle duplicates as needed
    else:
        return None

def df_to_msp(df, mspfile, outfile, newcol="RT", verbose=True):
    with open(mspfile, 'r') as f:
        lines = f.readlines()

    with open(outfile, 'w') as f:
        name = ""
        inchi = ""
        smiles = ""
        for i, line in enumerate(lines):
            # Print progress every 50000 lines
            if i % 50000 == 0:
                print(f"Processing line {i}... of {len(lines)}")
            if name == "":
                f.write(line)
                if "NAME:" in line:
                    name = line.split(":", 1)[1].strip()
            else:
                if "INCHIKEY:" in line:
                    inchi = line.split(":", 1)[1].strip()
                    f.write(line)
                    continue

                if "SMILES:" in line:
                    smiles = line.split(":", 1)[1].strip()
                    f.write(line)
                    continue

                if "RETENTIONTIME:" in line:
                    newrt = get_rt_for_name(df, name, inchi, smiles, newcol, verbose=verbose)
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

if __name__ == "__main__":
    # os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")
    # # posdf = msp_to_df("MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp")
    # # # Write to CSV
    # # posdf.to_csv("Pos_C1_Ref_FromMSP.csv", index=False)
    # #
    # # # Read neg file
    # # negdf = msp_to_df("MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp")
    # # # Write to CSV
    # # negdf.to_csv("Neg_C1_Ref_FromMSP.csv", index=False)
    #
    # # Test the df to msp function
    # posdf = pd.read_csv("Pos_C1_Ref_FromMSP.csv")
    # df_to_msp(posdf, "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp", "Test_Pos_C1.msp")
    #
    # exit()

    os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

    posfile = "PosIDsC1_Ref.csv"
    negfile = "NegIDsC1_Ref.csv"

    posdf = pd.read_csv(posfile)
    negdf = pd.read_csv(negfile)
    # Combine dataframes
    df = pd.concat([negdf, posdf])

    # Add Tail Length and Unsaturation columns if not present
    # if "Tail Length" not in df.columns:
    #     df = LipidFunctions.clean_df(df, column="Metabolite name")
    #     df = LipidFunctions.set_tails(df, columnname="Metabolite name")

    df = df_to_df(df)

    outfile = "C1_Ref_Training.csv"
    df.to_csv(outfile, index=False)

### Notes
"""
Used RetipFormatter to take output from MSDial and format it for Retip training.
Ran RetipPreprocessor to prepare the data for training
Did initial training for 30 min with S/N average as weight column and the top 10 features plus nAcid and nBase as the descriptor names.
Used high quality preset for training.
Ran RetipPredictor to predict RTs for the training set again. 
Used TailTrendsPlotter to analyze residuals and filter out entries with residuals above 0.25 min.
Reran RetipTrainer on the filtered data for 10 min to get the final model with high quality and no weight column.
Ran RetipLibProcessor to process MSP file into CSV for prediction and add descriptors
Ran RetipLibPredictor to predict RTs for the MSP entries in the CSV
Ran RetipLibtoMSP to write the predicted RTs back into a new MSP file.

"""

