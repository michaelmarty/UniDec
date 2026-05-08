from PublicScripts.Lipids.LipidFunctions import *
import os

topdir = r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs"
os.chdir(topdir)

# df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
df = pd.read_excel("Molecule Transition Results_stats.xlsx", sheet_name="Molecule Stats")

print(len(df))
mzrt_plot(df, write_file=True, class_col="Molecule List Name", rtcol="RT", mzcol="MZ")
