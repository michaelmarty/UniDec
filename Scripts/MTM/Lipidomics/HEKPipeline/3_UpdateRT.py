import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2"
os.chdir(topdir)

tlpath = r"Combined_TransitionList_Cleanest_4.csv"
rtfile = r"RT_Shifting.tsv"
output_path = tlpath.replace("_4.csv", "_5.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")
forcedrts = {"CL 70:3|CL 16:0_18:1_18:1_18:1": 15.6, "CL 70:3|CL 16:1_18:0_18:1_18:1": 15.6,
             "CL 72:3|CL 18:0_18:1_18:1_18:1": 15.6, "TG 64:0|TG 16:0_18:0_30:0": 18.5, "TG 68:1|TG 16:0_18:0_34:1": 18.8,
             "TG 68:2|TG 16:0_18:1_34:1": 18.7, "TG 68:2|TG 18:0_18:1_32:1": 18.5, "TG 70:2|TG 18:0_34:0_18:2": 18.6}
df = update_rt(df, cdf, mode="apex", min_window=0.5, max_window=2.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


