import os
import pandas as pd
import numpy as np
# from LipidFunctions import lipidlist_to_LC
from LipidFunctions import *


topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

spectrumkey = "Ref Spec" # or "MS/MS spectrum" if experimental spectra are used, if so, a cutoff is recommended
cutoff=0

infile = "NegIDsC1_Ref.csv"
outfile = "Neg_TransitionList.csv"

df = pd.read_csv(infile)
df = msdial_to_transitionlist(df, cutoff=cutoff, spectrumkey=spectrumkey)
df.to_csv(outfile, index=False)

# Repeate for Positive Mode
infile = "PosIDsC1_Ref.csv"
outfile = "Pos_TransitionList.csv"

posdf = pd.read_csv(infile)
posdf = msdial_to_transitionlist(posdf, cutoff=cutoff, spectrumkey=spectrumkey)
posdf.to_csv(outfile, index=False)

# Merge pos and neg
df = pd.concat([df, posdf])
# Write to csv
df.to_csv("Combined_TransitionList.csv", index=False)
# Print number of unique Molecule Names
print("Number of unique Molecule Names:", df["Molecule Name"].nunique())

# Needed: Precursor Ion Formula	Precursor Neutral Formula
# Needed: Product Mz	Product Charge	Product Ion Formula	Product Neutral Formula	Product Adduct