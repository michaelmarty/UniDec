import os
import pandas as pd
import numpy as np
from LipidFunctions import merge_pos_neg

if __name__ == "__main__":
    topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
    os.chdir(topdir)

    posfile = "PosIDsC1_Ref.csv"
    negfile = "NegIDsC1_Ref.csv"

    posdf = pd.read_csv(posfile)
    negdf = pd.read_csv(negfile)

    rttol = 0.1

    newposdf, newnegdf = merge_pos_neg(posdf, negdf, rttol)

    newposdf.to_csv("PosIDsC1_Ref.csv", index=False)
    newnegdf.to_csv("NegIDsC1_Ref.csv", index=False)

