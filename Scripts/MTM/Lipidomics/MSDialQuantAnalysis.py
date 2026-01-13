import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import *


topdir = r"Z:\Group Share\Annika\Stellar\FAM\Analysis"

set1 = ["20260108_AS_02_AcidifiedHEK_UntDDA", "20260108_AS_09_AcidifiedHEK_UntDDA", "20260108_AS_10_AcidifiedHEK_UntDDA"]
set2 = ["20260108_AS_08_CB2_disc_UntDDA", "20260108_AS_11_CB2_disc_UntDDA", "20260108_AS_12_CB2_disc_UntDDA"]
set2 = ["20260108_AS_08_CB2_disc_UntDDA", "20260108_AS_11_CB2_disc_UntDDA"]

valfile = "PeakValues_3_2026_01_13_09_54_37.csv"
masterfile = "PeakMaster_3_2026_01_13_09_54_37.csv"

vdf = pd.read_csv(f"{topdir}\\{valfile}")
mdf = pd.read_csv(f"{topdir}\\{masterfile}")
sdfpos = translate_to_skyline(vdf, mdf)

negvalfile = "PeakValues_2_2026_01_13_11_11_49.csv"
negmasterfile = "PeakMaster_2_2026_01_13_11_11_49.csv"
vdfneg = pd.read_csv(f"{topdir}\\{negvalfile}")
mdfneg = pd.read_csv(f"{topdir}\\{negmasterfile}")
sdfneg = translate_to_skyline(vdfneg, mdfneg)

sdf = pd.concat([sdfpos, sdfneg], ignore_index=True)

# Drop FA class
sdf = sdf[sdf['Molecule List Name'] != 'FA']

# Write to file
file = f"{topdir}\\SkylineInput.csv"
sdf.to_csv(file, index=False)

full_pipeline(file, set1, set2, mode="All", drop_IS=False, normalize_IS=False,
              normalize_TIC=True, paired=False, bh_correction=True, plot_results=True, heatplots=False)

plt.savefig(f"{topdir}\\Lipidomics_Comparison.pdf")