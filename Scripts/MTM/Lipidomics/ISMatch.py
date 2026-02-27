import os
import pandas as pd
from PublicScripts.Lipids.SkylineResultsPlotter import default_isdf as isdf
from PublicScripts.Lipids.LipidFunctions import gen_smiles_similarity_matrix
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Wxagg")
import numpy as np

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

def get_gaussian_similarity_matrix(mzlist1, mzlist2, sd=50.0, shift=7.0):
    matrix = np.zeros((len(mzlist1), len(mzlist2)))
    for i, mz1 in enumerate(mzlist1):
        for j, mz2 in enumerate(mzlist2):
            # Gaussian similarity based on difference in m/z
            diff = abs(mz1 - (mz2 - shift))
            matrix[i, j] = np.exp(-0.5 * (diff / sd)**2)
    return matrix


topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Shivam\IQX_Stellar_fixed"
os.chdir(topdir)

filepath = r"Combined_DDA_OutlierAnalysis_Results.xlsx"
df = pd.read_excel(filepath, sheet_name="Predicted Good")

# Isolate Dataset is Stellar
df = df[df["Dataset"] == "Stellar"]
# Isolate CL lipids
# df = df[df["Ontology"] == "PE"]
# Reset df index
df = df.reset_index(drop=True)

# Isolate PC class in isdf
# isdf = isdf[isdf["Class"] == "PC"]
# isdf = isdf.reset_index(drop=True)

slist = df["SMILES"].tolist()
islist = isdf["SMILES"].tolist()

matrix = gen_smiles_similarity_matrix(slist, islist, simtype="tanimoto", fingerprinttype="rdkit")
print(matrix.shape)
#
# mzlist = df["Average Mz"].tolist()
# rtlist = df["Average Rt(min)"].tolist()
#
# ismzlist = isdf["Mass"].tolist()
# isrtlist = isdf["RT"].tolist()
#
# mzmatrix = get_gaussian_similarity_matrix(mzlist, ismzlist, sd=100, shift=7)
# rtmatrix = get_gaussian_similarity_matrix(rtlist, isrtlist, sd=0.5, shift=0)
#
# matrix = matrix * mzmatrix

# Normalize each row of the matrix
row_sums = matrix.max(axis=1, keepdims=True)
# Avoid division by zero
row_sums[row_sums == 0] = 1
matrix = matrix / row_sums

# Add each row of the matrix into the dataframe
for i, row in df.iterrows():
    vector = matrix[i]
    print(row["Metabolite name"], isdf.at[vector.argmax(), "Lipid"])
    df.at[i, "Max Similarity"] = max(vector)
    df.at[i, "Best Match"] = isdf.at[vector.argmax(), "Lipid"]
    df.at[i, "Best Match Index"] = vector.argmax()
    df.at[i, "SimVector"] = str(vector)
    # plot_smiles(row["SMILES"], isdf.at[vector.argmax(), "SMILES"])

# print(df.to_string())
plt.figure(figsize=(30, 12))
plt.imshow(matrix.transpose(), aspect="auto", cmap="viridis")
plt.xlabel("Predicted Good Lipids")
plt.ylabel("ISDF Lipids")
# Only label x ticks for classes
unique_classes = df["Ontology"].unique()
class_indices = {cls: df[df["Ontology"] == cls].index[0] for cls in unique_classes}
plt.xticks(list(class_indices.values()), list(class_indices.keys()), rotation=90)
plt.yticks(range(len(isdf)), isdf["Lipid"], rotation=0)

plt.tight_layout()
plt.show()

# Save the updated dataframe to a new Excel file
# output_filepath = r"Combined_DDA_OutlierAnalysis_Results_with_Similarity.xlsx"
# df.to_excel(output_filepath, index=False)
