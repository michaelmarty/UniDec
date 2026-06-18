import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import tl_to_il, derplicate_il

topdir = r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\20260611"
os.chdir(topdir)

tlpath = r"HEK_TransList_Clean_CB2_2.csv"
outpath = r"CB2_IsoList_Clean.csv"
tldf = pd.read_csv(tlpath)

ildf = tl_to_il(tldf)
print(len(ildf))
print("Number of unique names:", ildf["Compound"].nunique())
ildf = derplicate_il(ildf)
print(len(ildf))

ildf.to_csv(outpath, index=False)

