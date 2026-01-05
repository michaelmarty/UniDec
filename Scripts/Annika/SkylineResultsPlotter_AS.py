import os
import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import full_pipeline

if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\EC\Test_TransResults.csv"

    set1 = ["20251216_AS_11_FT1_EC", "20251216_AS_13_FT2_EC", "20251216_AS_15_FT3_EC"]
    set2 = ["20251216_AS_05_D1_EC_2", "20251216_AS_07_D2_2", "20251216_AS_09_D3_EC_2"]

    full_pipeline(file, set1, set2, mode="Products", drop_IS=True, normalize_IS=False,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True)
