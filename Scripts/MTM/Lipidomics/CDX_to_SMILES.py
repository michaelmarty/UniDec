import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Wxagg")

def plot_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    img = Draw.MolToImage(m, size=(900, 900))
    plt.imshow(img)
    plt.axis("off")
    plt.show()

path = r"C:\Data\ChemDrawStds\700041.cdx"
rdkit_mol = Chem.MolsFromCDXMLFile(path)
rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
print(rdkit_smiles)
plot_smiles(rdkit_smiles)
exit()

os.chdir(r"C:\Data\ChemDrawStds")

# Read csv file with product numbers and names
csv_path = r"C:\Data\ChemDrawStds\IS_Information_with_SMILES.csv"
isdf = pd.read_csv(csv_path)

# Loop through smiles and plot
for index, row in isdf.iterrows():
    smiles = row["SMILES"]
    plot_smiles(smiles)

exit()


exe = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"

# List all files with .cdx extension
cdx_files = [f for f in os.listdir() if f.endswith(".cdx")]

outdict = {}

for f in cdx_files:
    number = f.split(".")[0]
    # If file exists, skip
    if not os.path.exists(f"{number}.smi"):
        call = f'"{exe}" -i cdx "{f}" -o smi -O "{number}.smi"'
        print(f"Processing {f} with command: {call}")
        subprocess.run(call, shell=True)

    # Open smi file and read the SMILES string
    with open(f"{number}.smi", "r") as smi_file:
        smiles = smi_file.read().strip()
        # If "[Na+]" in the SMILES, drop it, same with NH4
        if "[Na+]" in smiles:
            smiles = smiles.replace("[Na+]", "")
        if "[NH4+]" in smiles:
            smiles = smiles.replace("[NH4+]", "")
        smiles = smiles.strip()
        print(f"SMILES for {f}: {smiles}")
        outdict[number] = smiles
        plot_smiles(smiles)

# Create a dataframe from the output dictionary
df = pd.DataFrame(list(outdict.items()), columns=["Product Number", "SMILES"])

# Merge with this dataframe
ispath = r"Z:\Group Share\Annika\Stellar\IS_Information_2.csv"
isdf = pd.read_csv(ispath)

# Set product number to int in both dataframes for merging
isdf["Product Number"] = isdf["Product Number"].astype(int)
df["Product Number"] = df["Product Number"].astype(int)

merged_df = pd.merge(isdf, df, on="Product Number", how="left")

# Write out the merged dataframe to a new CSV file
merged_df.to_csv("IS_Information_with_SMILES.csv", index=False)

