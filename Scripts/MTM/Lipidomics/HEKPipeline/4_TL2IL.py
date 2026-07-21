import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import tl_to_il, derplicate_il

topdir = r"Z:\Group Share\Estevan\Lipid Stuffs"
os.chdir(topdir)

tlpath = r"HEK_TransList_Clean3.csv"
outpath = r"HEK_IsoList_Clean3.csv"
tldf = pd.read_csv(tlpath)


ildf = tl_to_il(tldf)
print(len(ildf))
print("Number of unique names:", ildf["Compound"].nunique())
ildf = derplicate_il(ildf)
print(len(ildf))

ildf.to_csv(outpath, index=False)

