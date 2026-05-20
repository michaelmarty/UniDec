from PublicScripts.Lipids.LipidFunctions import *
import os

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData"
os.chdir(topdir)

df = pd.read_csv("HEK_TransList_Clean.csv")

print(len(df))
print(len(df["Molecule Name"].unique()))
print(len(df["Molecule List Name"].unique()))
mzrt_plot(df, write_file=True, class_col="Molecule List Name", rtcol="Explicit Retention Time", mzcol="Precursor Mz")
