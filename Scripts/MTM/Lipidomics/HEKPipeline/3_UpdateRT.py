import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import update_rt

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Final_results"
os.chdir(topdir)

tlpath = r"Combined_TransitionList_Clean.csv"
rtfile = r"HEK.tsv"
output_path = tlpath.replace(".csv", "_RT.csv")

df = pd.read_csv(tlpath)

cdf = pd.read_csv(rtfile, sep="\t")
forcedrts = {"CL 68:2|CL 16:0_18:1_18:0_16:1": 15.42}
df = update_rt(df, cdf, mode="apex", min_window=3.0, max_window=5.0, nstds=8, forcedrts=forcedrts)

# Write to output
df.to_csv(output_path, index=False)


