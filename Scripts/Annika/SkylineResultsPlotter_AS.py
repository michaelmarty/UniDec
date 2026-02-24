import os
import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import cv_pipeline

if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Extr_US_MoleculeTransResults.csv"

    # set1 = ["0.5% FT1", "0.5% FT2", "0.5% FT3"]
    set2 = ["08", "09", "10"]

    cv_pipeline(file, results_set=set2, mode="Products", drop_IS=False, normalize_IS=True, normalize_TIC=True,)

    exit()