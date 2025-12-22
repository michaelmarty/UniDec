import os
import pandas as pd
from PublicScripts.Lipids.LipidFragPrediction import assign_df

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

negfile = "PosIDsC1_Ref.csv"

negdf = pd.read_csv(negfile)
# Filter negdf to only things where Metabolite name is PG 18:1_18:1
# negdf = negdf[negdf["Metabolite name"].str.contains("PG 18:1_18:1")]
outdf = assign_df(negdf)

outdf.to_csv(negfile[:-4] + "_TailCheck.csv", index=False)
