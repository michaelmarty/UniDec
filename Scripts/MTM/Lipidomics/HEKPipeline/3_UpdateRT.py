import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2"
os.chdir(topdir)

tlpath = r"AllDatasets_TransList_Clean3.csv"
rtfile = r"RT_Shifts.tsv"
output_path = tlpath.replace("_Clean3.csv", "_Clean4.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")
forcedrts = {"CL 70:2|CL 16:0_18:1_18:0_18:1": 14.2,"TG 64:0|TG 16:0_18:0_30:0": 18.5, "TG 68:1|TG 16:0_18:0_34:1": 18.8, "CL 66:1|CL 16:0_16:0_16:0_18:1": 14.1,
             "TG 68:2|TG 16:0_18:1_34:1": 18.7, "TG 68:2|TG 18:0_18:1_32:1": 18.5, "TG 70:2|TG 18:0_34:0_18:2": 18.6, "CL 70:2|CL 16:0_18:1_18:0_18:1": 14.3,
             "CL 72:3|CL 18:0_18:1_18:1_18:1": 14.3, "CL 70:3|CL 16:1_18:0_18:1_18:1": 14.2, "CL 66:2|CL 16:0_16:1_16:0_18:1": 14.0, "CL 68:2|CL 16:0_18:1_16:0_18:1": 14.0,
             "CL 70:4|CL 16:1_18:1_18:1_18:1": 13.9, "Cer 42:1;2|Cer 18:1;2/24:0": 13.6, "Cer 32:1;2|Cer 18:1;2/14:0": 7.8, "Cer 34:0;2|Cer 18:0;2/16:0_A": 9.3, "Cer 37:3;2|Cer 18:1;2/19:2": 9.8,
             "Cer 40:0;2|Cer 18:0;2/22:0_A": 13.4, "Cer 40:1;2|Cer 18:1;2/22:0_A": 13.4, "Cer 41:2;2|Cer 19:1;2/22:1": 12.3, "Cer 42:1;2|Cer 18:1;2/24:0_A": 13.7, "Cer 42:1;2|Cer 18:1;2/24:0_B": 13.7,
             "Cer 42:2;2|Cer 18:1;2/24:1_A": 13.2, "Cer 48:2;2|Cer 18:1;2/30:1": 14, "Cer 42:0;3|Cer 18:0;3/24:0": 13.5, "Cer 44:3;4|Cer 31:2;3/13:1": 13, "DG 35:2|DG 17:1_18:1": 11,
                "FAHFA 34:2;O|FAHFA 16:0/18:1;O": 8.4, "FAHFA 35:8;O|FAHFA 18:5/17:2;O": 8.9, "LPC 14:1": 2.2, "LPC 22:0": 5.8, "LPC 24:0": 6.6, "LPC 26:1": 6.6, "LPC O_22:2": 5.3,
             "LPG 14:0": 2.2, "PC 30:1|PC 14:0_16:1": 6.8, "PC 33:0|PC 16:0_17:0_B": 8.9, "PC 33:0|PC 16:0_17:0_A": 8.9, "PC 38:2|PC 18:1_20:1": 10.1, "PC O_38:2_A": 9.4, "PC O_38:2_B": 9.4,
             "PC O_40:1|PC O_22:0_18:1": 13.5, "PE 32:2|PE 16:1_16:1": 7.1, "PE 38:3|PE 18:0_20:3": 10.1, "PE 40:3|PE 18:0_22:3": 11.2, "PG 38:4|PG 18:1_20:3": 6.5, "PI 32:1|PI 16:0_16:1": 6.4,
             "PI 32:2|PI 14:1_18:1": 7, "PI 36:3|PI 16:0_20:3": 7, "PI 36:3|PI 18:1_18:2": 6.7, "PI 38:4|PI 18:1_20:3": 7, "PI O_32:1|PI O_16:0_16:1": 6.9, "PI O_32:2|PI O_18:1_14:1": 6.3,
             "PI O_34:1|PI O_16:0_18:1": 7.5, "PI O_34:2": 7.1, "PI O_36:1|PI O_18:0_18:1": 8.5, "PI O_38:2|PI O_18:1_20:1": 8.5, "PS 30:1|PS 14:0_16:1": 5.7, "PS 36:3|PS 18:1_18:2": 6.8,
             "PS 38:1|PS 20:0_18:1": 9.2, "PS 38:3|PS 18:1_20:2": 7.5, "PS 42:1|PS 24:0_18:1": 14.9, "SM 40:1;2|SM 18:1;2/22:0": 11.5, "TG 68:2|TG 16:0_18:1_34:1": 16.7

             }
df = update_rt(df, cdf, mode="apex", min_window=0.5, max_window=2.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


