import pandas as pd
import os
from io import StringIO
from PublicScripts.Lipids.LipidFragPrediction import *

# Define CSV

lipid_data = """
Class	Adduct	Required Frags
Cer	[M+H]	FA1, FA2
Cer	[M-H]	FA1, FA2
Cer	[M+CH3COO]	FA1, FA2
CL 	[M-H]	FA1, FA2, FA3, FA4
CoQ	[M+H]	HG
DG	[M+NH4]	FA1, FA2
FAHFA	[M-H]	FA1
GM3	[M-H]	HG
HexCer	[M+H]	HG, FA1
HexCer	[M-H]	HG, FA1
HexCer	[M+CH3COO]	HG, FA1
LPA	[M-H]	HG
LPC	[M+H]	HG
LPC O	[M+CH3COO]	HG
LPC	[M+CH3COO]	HG, FA1
LPC O	[M+H]	HG
LPE	[M-H]	HG, FA1
LPG	[M-H]	HG, FA1
LPI	[M-H]	HG, FA1
LPS	[M-H]	HG, FA1
PA	[M-H]	HG, FA1, FA2
PC	[M+H]	HG(184)
PC O	[M+CH3COO]	HG, FA1, FA2
PC	[M+CH3COO]	HG, FA1, FA2
PC O	[M+H]	HG(184)
PE	[M+H]	HG, FA1, FA2
PE	[M-H]	HG, FA1, FA2
PE O	[M+H]	HG, FA1, FA2
PE O	[M-H]	HG, FA1, FA2
PG	[M-H]	HG, FA1, FA2
PI	[M-H]	HG, FA1, FA2
PS	[M-H]	HG, FA1, FA2
SM	[M+CH3COO]	HG
SM	[M+H]	HG(184)
TG	[M+NH4]	FA1, FA2, FA3

"""
# Read the CSV data into a pandas DataFrame

lipid_df = pd.read_csv(StringIO(lipid_data), sep="\t")



if __name__ == "__main__":
    filepath = r"Z:\Group Share\Annika\DDA_Compare\mzmine\Combined_DDA_Full.xlsx"
    os.chdir(os.path.dirname(filepath))
    df = pd.read_excel(filepath)

    # Filter to Ontology == PC
    df = df[df["Ontology"] == "PC"]
    print("Filtered PC DataFrame:")

    df = assign_df_fragments()

    # Write df to test output
    df.to_excel("PC_Fragment_Assignment.xlsx", index=False)