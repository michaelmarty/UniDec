import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import tl_to_il, derplicate_il

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2"
os.chdir(topdir)

tlpath = r"AllDatasets_TransList_Clean4.csv"
outpath = r"AllDatasets_IsoList_Clean4.csv"

tldf = pd.read_csv(tlpath)

ildf = tl_to_il(tldf)
print(len(ildf))
print("Number of unique names:", ildf["Compound"].nunique())
ildf = derplicate_il(ildf)
print(len(ildf))

ildf.to_csv(outpath, index=False)

