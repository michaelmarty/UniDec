import os
import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import full_pipeline


if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\20251216_MoreLipidDiscs_FTs\DDM Comparisons\0pt5DDM_Discsvs_FTs.csv"

    set1 = ["0.5% FT1", "0.5% FT2", "0.5% FT3"]
    set2 = ["0.5% D1", "0.5% D2", "0.5% D3"]

    full_pipeline(file, set1, set2, mode="Products", drop_IS=True, normalize_IS=True,
                  normalize_TIC=False, paired=True, bh_correction=True, plot_results=True)
