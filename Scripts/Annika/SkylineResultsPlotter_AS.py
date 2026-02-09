import os
import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import full_pipeline


if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\EC\0.5% DDM\Molecule Ratio Results.csv"

    set1 = ["0.5% FT1", "0.5% FT2", "0.5% FT3"]
    set2 = ["0.5% D1", "0.5% D2", "0.5% D3"]

    full_pipeline(file, set1, set2, mode="Products", drop_IS=False, normalize_IS=False,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True)
