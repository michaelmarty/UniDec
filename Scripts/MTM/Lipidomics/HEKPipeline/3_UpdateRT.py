import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Estevan\Lipid Stuffs"
os.chdir(topdir)

tlpath = r"HEK_TransList_Clean2.csv"
rtfile = r"RT_Shifts.tsv"
output_path = tlpath.replace("_Clean2.csv", "_Clean3.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")


# forcedrts = {"Cer 44:1;2|Cer 18:1;2/26:0": 11.3, "Cer 50:2;2|Cer 18:1;2/32:1": 12, "Cer 50:4;2|Cer 18:1;2/32:3": 12, "PC 30:1|PC 14:0_16:1": 5.7,
#              "PC 30:2|PC 14:1_16:1": 5.0,  "PC 33:0|PC 16:0_17:0_B": 8.3, "PC 34:0|PC 16:0_18:0": 7.7, "PC 40:0|PC 16:0_24:0": 9.8, "PE 40:4|PE 18:1_22:3": 7.5,
#              "PG 30:1|PG 14:0_16:1": 4.5, "IS PI(+[2]H5) 16:1_17:0": 5.3, "IS PI(+[2]H5) 17:0_18:1": 5.9, "PI 32:2|PI 14:1_18:1": 4.5, "PS 34:2|PS 16:1_18:1": 5.2,
#              "PS 38:3|PS 18:1_20:2": 5.9, "SM 31:3;2": 6.7, "SM 34:0;2": 6.2, "IS TG(+[2]H5) 13:0_14:0_14:0": 13, "IS TG(+[2]H5) 14:0_14:0_15:1": 13, "IS TG(+[2]H5) 14:0_14:0_17:1": 13,
#              "IS TG(+[2]H5) 15:1_16:0_16:0": 13, "PC O_34:0|PC O_18:0_16:0": 8.3, "PE O_36:0|PE O_16:0_20:0": 9.5, "PE O_38:2|PE O_18:1_20:1": 8.5, "PE 22:1_18:1": 9.5,
#
# }
forcedrts = {"Cer 44:1;2|Cer 18:1;2/26:0": 12.0, "Cer 44:1;2|Cer 20:1;2/24:0": 12.0, "HexCer 34:1;3|HexCer 18:1;2/16:0;O": 4.8, "HexCer 42:1;2|HexCer 18:1;2/24:0": 10.3,
             "PE 40:3|PE 18:0_22:3": 8.4, "PE 40:4|PE 18:1_22:3": 7.5, "IS PI(+[2]H5) 16:1_17:0": 5.3, "IS PI(+[2]H5) 17:0_18:1": 5.8, "PS 38:3|PS 18:1_20:2": 5.9,
             "PS 34:2|PS 16:1_18:1": 5.3, "PS 38:2|PS 18:1_20:1_A": 6.3, "PC O_34:0|PC O_18:0_16:0": 8.4, "LPI 18:1": 2.1, "LPI 20:1": 2.4, "IS PC(+[2]H5) 17:0_22:4": 6.9, "PC 30:1|PC 14:0_16:1": 5.5,
             "PC 33:0|PC 16:0_17:0_A": 5.8, "PC 34:0|PC 16:0_18:0": 7.5, "PC 40:0|PC 16:0_24:0": 9.4, "PC 42:0|PC 16:0_26:0": 12.2, "PE 35:1|PE 17:0_18:1": 7.3, "PE 35:2|PE 17:1_18:1": 6.6,
             "PE 40:4|PE 18:1_22:3": 7.3, "PG 30:1|PG 14:0_16:1": 4.4, "IS PI(+[2]H5) 16:1_17:0": 5.1, "IS PI(+[2]H5) 17:0_18:1": 5.7, "PI 32:2|PI 14:1_18:1": 4.4, "PS 38:3|PS 18:1_20:2": 5.7,
             "IS TG(+[2]H5) 13:0_14:0_14:0": 13}



forcedrts = {
             }
df = update_rt(df, cdf, mode="apex", min_window=1, max_window=2.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


