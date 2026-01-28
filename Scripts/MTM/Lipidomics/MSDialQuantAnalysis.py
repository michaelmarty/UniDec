import pandas as pd
import numpy as np
from PublicScripts.Lipids.SkylineResultsPlotter import *

def norm_pos_neg(posdf, negdf, common=True):
    if common:
        # Find common species in posdf and negdf
        common_species = set(posdf['Molecule']).intersection(set(negdf['Molecule']))
        ratios = []
        for species in common_species:
            pos_area = posdf[posdf['Molecule'] == species]['Area'].values[0]
            neg_area = negdf[negdf['Molecule'] == species]['Area'].values[0]
            if neg_area > 0:
                ratio = pos_area / neg_area
                ratios.append(ratio)
        if len(ratios) == 0:
            print("No common species with non-zero area found for normalization.")
            exit()
        avg_ratio = np.median(ratios)
        # Normalize negdf areas by avg_ratio
        negdf['Area'] = negdf['Area'] * avg_ratio
        print("Normalized negative mode by average ratio of common species:", avg_ratio)
    else:
        # Normalize positive and negative dataframes to their own max
        posmax = posdf['Area'].max()
        negmax = negdf['Area'].max()
        posdf['Area'] = posdf['Area'] / posmax
        negdf['Area'] = negdf['Area'] / negmax
    return posdf, negdf


topdir = r"Z:\Group Share\Annika\Stellar\FAM\Analysis"

set1 = ["20260108_AS_02_AcidifiedHEK_UntDDA", "20260108_AS_09_AcidifiedHEK_UntDDA", "20260108_AS_10_AcidifiedHEK_UntDDA"]
set2 = ["20260108_AS_08_CB2_disc_UntDDA", "20260108_AS_11_CB2_disc_UntDDA", "20260108_AS_12_CB2_disc_UntDDA"]
set2 = ["20260108_AS_08_CB2_disc_UntDDA", "20260108_AS_11_CB2_disc_UntDDA"]

valfile = "PeakValues_3_2026_01_13_09_54_37.csv"
masterfile = "PeakMaster_3_2026_01_13_09_54_37.csv"
masterfile = "PeakMaster_1_2026_01_13_13_49_35.csv"
valfile = "PeakValues_1_2026_01_13_13_49_35.csv"

vdf = pd.read_csv(f"{topdir}\\{valfile}")
mdf = pd.read_csv(f"{topdir}\\{masterfile}")
sdfpos = translate_to_skyline(vdf, mdf)

negvalfile = "PeakValues_2_2026_01_13_11_11_49.csv"
negmasterfile = "PeakMaster_2_2026_01_13_11_11_49.csv"
vdfneg = pd.read_csv(f"{topdir}\\{negvalfile}")
mdfneg = pd.read_csv(f"{topdir}\\{negmasterfile}")
sdfneg = translate_to_skyline(vdfneg, mdfneg)

sdfpos, sdfneg = norm_pos_neg(sdfpos, sdfneg)
sdf = pd.concat([sdfpos, sdfneg], ignore_index=True)
# sdf = sdfpos
# Sorty by molecule class and name
sdf = sdf.sort_values(by=['Molecule List Name', 'Molecule'])

# Drop FA class
sdf = sdf[sdf['Molecule List Name'] != 'FA']
sdf = sdf[sdf['Molecule List Name'] != 'TG']

# Write to file
file = f"{topdir}\\SkylineInput.csv"
sdf.to_csv(file, index=False)

figsave = f"{topdir}\\Lipidomics_Comparison.pdf"

full_pipeline(file, set1, set2, mode="All", drop_IS=False, normalize_IS=False,
              normalize_TIC=True, paired=False, bh_correction=True, plot_results=True, heatplots=False, figsavepath=figsave)
