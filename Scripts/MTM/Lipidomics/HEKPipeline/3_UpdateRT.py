import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\RT Shifting"
os.chdir(topdir)

tlpath = r"Combined_DDA_TransitionLists.csv"
rtfile = r"RT_Shifts.tsv"
output_path = tlpath.replace("_TransitionLists.csv", "_TransitionLists_Updated.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")




forcedrts = {"Cer 34:0;O2|Cer 18:0;O2/16:0": 7.6, "Cer 42:0;O2|Cer 18:0;O2/24:0": 12.8, "Cer 32:1;O2|Cer 18:1;O2/14:0": 7, "Cer 34:1;O2|Cer 18:1;O2/16:0": 7.4, "Cer 34:2;O2|Cer 18:1;O2/16:1": 7,
             "Cer 36:1;O2|Cer 18:1;O2/18:0": 8.2, "Cer 42:1;O2|Cer 18:1;O2/24:0": 12.3, "Cer 42:2;O2|Cer 18:1;O2/24:1": 10.5, "Cer 44:2;O2|Cer 18:1;O2/26:1": 12.2, "Cer 48:2;O2|Cer 34:2;O2/14:0": 13.7,
             "Cer 52:3;O2|Cer 18:1;O2/34:2": 13.9, "LPE P-20:0": 4.5, "PC O-34:0|PC O-18:0_16:0": 8.7, "PC O-36:0|PC O-20:0_16:0": 10.0, "PC O-38:0|PC O-22:0_16:0": 11.4, "PE O-38:2|PE O-18:1_20:1": 8.9,
             "PE O-42:1|PE O-24:0_18:1": 13.1, "LPE-N (FA)29:1|LPE-N (FA 13:1)16:0": 5.6, "PS 34:3|PS 16:1_18:2": 4.9, "PS 38:1|PS 20:0_18:1": 7.0, "PS 36:3|PS 18:1_18:2": 5.5,
             }

forcedrts = {}
df = update_rt(df, cdf, mode="apex", min_window=1, max_window=2.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


