import os
import pandas as pd
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA_RT_Corrected.csv")

df = full_name_assignment(df, rttol=0.2, rtcol="Average Rt(min)")

df["Ontology"] = df["Ontology"].replace({
    # Combine these as Cer
    "Cer_ADS": "Cer",
    "Cer_AP": "Cer",
    "Cer_AS": "Cer",
    "Cer_EODS": "Cer",
    "Cer_EOS": "Cer",
    "Cer_NDS": "Cer",
    "Cer_NP": "Cer",
    "Cer_NS": "Cer",
    "CerP": "Cer",

    # Combine these as HexCer
    "AHexCer": "HexCer",
    "SHexCer": "HexCer",
    "HexCer_HDS": "HexCer",
    "HexCer_HS": "HexCer",
    "HexCer_NDS": "HexCer",
    "HexCer_NS": "HexCer",
})

df = get_number_of_datasets(df, namecol="Metabolite name")
write_dfs_to_excel(df)
