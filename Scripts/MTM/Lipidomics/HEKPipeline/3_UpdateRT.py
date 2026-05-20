import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData"
os.chdir(topdir)

tlpath = r"CleanedAdducts_AllDatasets_TransList_MM_Clean4.csv"
rtfile = r"RT_Shifts.tsv"
output_path = tlpath.replace("_Clean4.csv", "_Clean5.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")


forcedrts = {"Cer 34:0;2|Cer 18:0;2/16:0_A": 6.1, "Cer 37:3;2|Cer 18:1;2/19:2": 6.1,"HexCer 34:3;2|HexCer 18:2;2/16:1": 4.7, "LPC 14:1": 1.1, "LPE 17:1": 1.5, "LPI 22:2": 1.4, "PC 30:1|PC 14:0_16:1": 4.7,
             "PC 35:2|PC 17:1_18:1": 6.5, "PC 37:1": 6.0, "PC 37:2|PC 18:1_19:1": 6.3, "PC 38:2|PC 18:1_20:1": 6.7, "PC 38:3|PC 18:1_20:2": 6.2, "PC 40:0|PC 16:0_24:0": 8.0, "PC 48:1|PC 16:0_32:1": 13.1,
             "PE 38:0|PE 18:0_20:0": 7.5, "PE 38:4|PE 18:1_20:3": 6.3, "PE 40:3|PE 18:0_22:3": 7.1, "PE 42:2|PE 18:1_24:1": 8.3, "PG 30:1|PG 14:0_16:1": 3.3, "PG 32:1|PG 14:0_18:1": 4.0, "PG 34:1|PG 18:0_16:1": 4.6,
             "PG 34:2|PG 16:1_18:1": 4.1, "PI 32:1|PI 16:0_16:1": 3.9, "PI 32:0|PI 16:0_16:0": 4.4, "PI 32:2|PI 14:1_18:1": 2.8, "PI 34:0|PI 16:0_18:0": 4.5, "PI 34:2|PI 16:1_18:1": 4.0, "PI 36:2|PI 16:0_20:2": 4.5,
             "PI 40:2|PI 18:1_22:1": 5.6, "PS 34:2|PS 16:1_18:1": 4.1, "TG 49:2|TG 15:0_16:1_18:1": 13.3

}
forcedrts = {"Cer 41:2;2|Cer 19:1;2/22:1": 9.0, "Cer 44:1;2|Cer 18:1;2/26:0":11.6 , "LPI 20:1": 2.4, "LPS 20:0": 2, "PA 34:3|PA 16:1_18:2": 5.3, "PA 36:0|PA 18:0_18:0": 6.8,
             "PA 36:1|PA 18:0_18:1": 6.7, "PC 30:1|PC 14:0_16:1": 5.7, "PC 30:2|PC 14:1_16:1": 5.2, "PC 31:0|PC 15:0_16:0": 6.5, "PC 32:0|PC 16:0_16:0": 6.9, "PC 40:0|PC 16:0_24:0": 9.7,
             "PC 40:0|PC 16:0_24:0": 6.8, "PE 38:4|PE 18:1_20:3": 7.0, "PE 40:4|PE 18:1_22:3": 7.6, "PE 44:2|PE 18:1_26:1": 11.4, "PG 34:2|PG 16:1_18:1": 5.3,
             "PS 34:2|PS 16:1_18:1": 5.3, "PS 38:1|PS 20:0_18:1": 7.1, "IS SM(+[2]H9) 18:1;2/20:1": 7.0, "IS SM(+[2]H9) 18:1;2/22:1": 7.7, "SM 34:0;2": 6.5, "IS TG(+[2]H5) 13:0_14:0_14:0": 14,
             "IS TG(+[2]H5) 14:0_14:0_15:1": 14, "IS TG(+[2]H5) 14:0_14:0_17:1": 14, "IS TG(+[2]H5) 15:1_16:0_16:0": 14, "SHexCer 42:1;2|SHexCer 17:0;2/25:1": 7.0, "SHexCer 42:2;2|SHexCer 18:1;2/24:1": 6.4

             }
df = update_rt(df, cdf, mode="apex", min_window=0.5, max_window=2.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


