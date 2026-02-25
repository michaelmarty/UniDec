import os
import pandas as pd
from PublicScripts.Lipids.SkylineResultsPlotter import default_isdf as isdf
from PublicScripts.Lipids.LipidFunctions import gen_smiles_similarity_matrix
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Wxagg")

def plot_smiles(smiles1, smiles2):
    m = Chem.MolFromSmiles(smiles1)
    AllChem.Compute2DCoords(m)
    Chem.AddHs(m)
    img = Draw.MolToImage(m, size=(900, 900))

    m2 = Chem.MolFromSmiles(smiles2)
    AllChem.Compute2DCoords(m2)
    Chem.AddHs(m2)
    # Set isotopes to 0
    for atom in m2.GetAtoms():
        atom.SetIsotope(0)
    img2 = Draw.MolToImage(m2, size=(900, 900))
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(img)
    plt.axis("off")
    plt.subplot(1, 2, 2)
    plt.imshow(img2)
    plt.axis("off")
    plt.show()

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Shivam\IQX_Stellar_fixed"
os.chdir(topdir)

filepath = r"Combined_DDA_OutlierAnalysis_Results.xlsx"
df = pd.read_excel(filepath, sheet_name="Predicted Good")

# Isolate Dataset is Stellar
df = df[df["Dataset"] == "Stellar"]
# Isolate CL lipids
df = df[df["Ontology"] == "PC"]
# Reset df index
df = df.reset_index(drop=True)

# Isolate PC class in isdf
isdf = isdf[isdf["Class"] == "PC"]
isdf = isdf.reset_index(drop=True)

slist = df["SMILES"].tolist()
islist = isdf["SMILES"].tolist()
# In islist, replace * with C
islist = [s.replace("*", "C") for s in islist]

matrix = gen_smiles_similarity_matrix(slist, islist, simtype="tanimoto", fingerprinttype="rdkit")
print(matrix.shape)

# Add each row of the matrix into the dataframe
for i, row in df.iterrows():
    vector = matrix[i]
    print(vector)
    df.at[i, "Max Similarity"] = max(vector)
    df.at[i, "Best Match"] = isdf.at[vector.argmax(), "Lipid"]
    df.at[i, "Best Match Index"] = vector.argmax()
    df.at[i, "SimVector"] = str(vector)
    plot_smiles(row["SMILES"], isdf.at[vector.argmax(), "SMILES"])

print(df.to_string())

# Save the updated dataframe to a new Excel file
# output_filepath = r"Combined_DDA_OutlierAnalysis_Results_with_Similarity.xlsx"
# df.to_excel(output_filepath, index=False)
