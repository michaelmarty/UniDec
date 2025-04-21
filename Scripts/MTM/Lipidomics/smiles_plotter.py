from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Wxagg")

m = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCCN)OC(=O)CCCCCCCC=CCCCCCCCC")
#AllChem.Compute2DCoords(m)

fig = Draw.MolToMPL(m, fitImage=True)
plt.show()

