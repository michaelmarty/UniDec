import os
import pandas as pd
import numpy as np
import styleframe as sf
from unidec.modules.hramtools import calc_fromula_from_mass, reallysmallsearch

mass_c = 12.000000
mass_h = 1.00784
mass_o = 15.99491
mass_n = 14.00307
mass_p = 30.97376
mass_ch2 = mass_c + (mass_h * 2)
mass_h2o = (mass_h * 2) + mass_o
mass_glycerol = (mass_c * 3) + (mass_h * 8) + (mass_o * 3)
mass_phosphate = mass_p + (mass_o * 4) + mass_h
mass_formate = (mass_h + mass_c + (mass_o * 2))  # HCOO-
mass_acetate = (mass_c * 2) + (mass_h * 3) + (mass_o * 2)  # CH3COO-
mass_ch3 = (mass_c * 1) + (mass_h * 3)  # CH3
mass_d = 2.014 - mass_h  # Deuterium mass difference
# Sphingolipid fragment C2H7NO
mass_slfrag = (mass_c * 2) + (mass_h * 7) + mass_n + mass_o
# CER ref: https://link.springer.com/article/10.1007/s00216-013-7601-y

hgbasedict = {"PG": [228.04], "PI": [316.05], "SQDG": [282.04], "PS": [87.032, 223.0022], "OxPS": [87.034, 223.0022, -36.0225],
              "PC": [281.102], "OxPG": [228.04, -36.0225], "MGDG": [218.077], "EtherPG": [228.04 - mass_o, 194.03206, 210.02856],
              "EtherPI": [316.05 - mass_o, 316.05, 162.05], "EtherPS": [223.02716, 87.032 - mass_o, 223.0022 - mass_o, 87.032],
              "DGGA": [250.07], "Ac2PIM1": [460.1, 136.00], "DLCL": [210.03, 364.032],
              "AHexCer": [208.052, 222.08, 2.96706],
              "Hex2Cer": [338.0849, 176.0349, 294.05706, 204.02706, 308.07706, 218.04706],
              "Hex3Cer": [500.1349, 338.0849, 176.0349, 366.07706, 380.09706], "HexCer": [176.0349],
              "HexCer_HS": [193.04706, 207.06706, 12.9895], "HexCer_HDS": [193.04706, 207.06706, 12.9895],
              "HexCer_EOS": [193.04706, 207.06706, 12.9895, 191.0443, 158.0197, 176.0298, 204.02196, 190.05066],
              "Cer_AS": [28.97706, 103.02366, 117.04366, 89.04366, 57.02366, 42.98975],
              "Cer_NP": [31.04295, -59.03784, -61.04784, -43.03784, 42.99216, 28.98216],
              "Cer_AP": [89.04366, 135.05366, 71.0315],
              "Cer_BDS": [9.9848, 34.00206, 11.99994, -44.02434, -12.00434],
              "Cer_BS": [93.04246, 46.00116, -14.01724, 33.02276, -44.03154],
              "Cer_EBDS": [73.0528],

              "HexCer_AP": [71.0315, 212.05066, 182.04066, 237.08066, 283.09066, 149.02216, 179.03216, 204.06216,
                            250.07216, 190.05216, 144.04216, 119.01216, 89.00216, 236.05216],
              "ASM": [4.03, 4.04, 6.04, 6.05, 8.06, 8.07, 10.08, 12.09, 12.1, 62.03, 62.04, 64.05,
                      64.06, 66.06, 66.07, 68.08, 68.09, 70.1, 72.11, 72.12, 76.05, 76.06, 78.06, 78.07,
                      80.08, 80.09, 82.1, 84.11, 84.12, 86.13, 86.14], "Ac2PIM1": [324.10226, 460.1001],
              "EtherSMGDG": [298.04069], "EtherPE": [197.04, 179.03016],
              "EtherMGDG": [236.08660, 264.083060, 278.10016],
              "EtherDGDG": [426.13876, 380.1323, 398.1423],
              "EtherOxPC": [-36.0225],
              "EtherOxPE": [-36.0225, 197.0464, 179.03856],
              "OxPE":[-36.0225],
              "OxPC":[-36.0225],
              "OxPI":[-36.0225],
              "LNAPE": [179.0326, 197.0483], "LNAPS": [69.02086],
              "NAGly": [43.99206], "NAGlySer":[61.99796, 30.00796],
              "PI_Cer": [-3.04804],
              "PE_Cer": [-3.04804], "Cer_EOS": [39.997], "Cer_EODS": [39.997],
              "SL":[79.9594]}

ignore_fragments = {}
ignore_nls = {"Hex2Cer": [46.01, 46.0], "Hex3Cer": [46.01], "HexCer": [46.01], "AHexCer": [46.01], "SM":[60.02], "SL":[111.98]}
acelated_classes = ["AHexCer", "ASM"]
adduct_nls = {"[M+CH3COO]-": [60.02], "[M+HCOO]-": [46.01]}


def predict_tail_fragments(tail, hgfrags=[], hgbases=[], hgnls=[], adduct="[M-H]-", classname=""):
    # If hgfrags is a string, split by comma
    if isinstance(hgfrags, str):
        hgfrags = hgfrags.split(",")
    if not isinstance(hgfrags, list):
        hgfrags = [hgfrags]

    if isinstance(hgnls, str):
        hgnls = hgnls.split(",")
    if not isinstance(hgnls, list):
        hgnls = [hgnls]

    if isinstance(hgbases, str):
        hgbases = hgbases.split(",")
    if not isinstance(hgbases, list):
        hgbases = [hgbases]

    fragments = {}
    nls = {}
    if pd.isna(tail):
        return fragments, nls
    if ":" in tail:
        parts = tail.split(":")
        numox = 0
        nhydrox = 0
        nmethyl = 0
        ndeut = 0
        nbound = False
        obound = False
        if "methyl" in parts[1]:
            parts[1] = parts[1].split("(methyl")[0]
            nmethyl = 1

        if "(d7)" in parts[1]:
            parts[1] = parts[1].replace("(d7)", "")
            ndeut = 7
        if "(d9)" in parts[1]:
            parts[1] = parts[1].replace("(d9)", "")
            ndeut = 9

        if "O-" in parts[0]:
            parts[0] = parts[0].replace("O-", "")
            obound = True

        if "N-" in parts[0]:
            parts[0] = parts[0].replace("N-", "")
            nbound = True

        if ";" in parts[1]:
            parts[1] = parts[1].split(";")[0]
            oxstring = tail.split(";")[1]
            if "(FA" in oxstring:
                oxstring = oxstring.replace("(FA", "")
            if "OH)" in oxstring:
                if "2OH" in oxstring:
                    nhydrox = 2
                elif "3OH" in oxstring:
                    nhydrox = 3
                elif "OH" in oxstring:
                    nhydrox = 1
                else:
                    print("Error parsing hydroxyls in tail:", tail)
                    raise ValueError
                oxstring = ""

            # Drop O in oxstring
            if oxstring == "O":
                numox = 1
            elif oxstring == "2O":
                numox = 2
            elif oxstring == "3O":
                numox = 3
            elif oxstring == "":
                numox = 0
            else:
                try:
                    numox = int(oxstring.replace("O", ""))
                except ValueError:
                    print("Error parsing oxidation in tail:", tail)
                    print(oxstring)
                    numox = 0

        try:
            num_carbons = int(parts[0])
            num_double_bonds = int(parts[1])
        except ValueError:
            print("Error parsing tail:", tail)
            return fragments, nls

        acylmass = num_carbons * mass_ch2 - mass_h + 2 * mass_o - (num_double_bonds * 2 * mass_h)

        # Predict common fragments
        mass = acylmass
        if numox > 0 and "Ox" not in classname:
            mass += mass_n + 3 * mass_h
            if numox == 1:
                mass -= mass_o
            if numox == 3:
                mass += mass_o
        elif numox > 0 and "Ox" in classname:
            mass += numox * mass_o
        # mass += nhydrox * (-mass_h)  # Add hydroxyls
        mass += nmethyl * (mass_c + 2 * mass_h)  # Add methyl group

        # If O bound, remove O and add 2 H
        if obound and classname not in acelated_classes:
            mass += - mass_o + 2 * mass_h

        # print(numox, nhydrox, nmethyl, ndeut, num_carbons, num_double_bonds)

        # Compute fragment and NL values without rounding; round them all at the end
        fragments["FA(+O)"] = mass

        # For Ceramides
        fragments["FA(-2H)"] = mass - mass_h2o

        fragments["FA(+NH)"] = mass - mass_o + mass_n + mass_h
        fragments["FA(+H2O)"] = mass - mass_o + mass_h2o
        fragments["FA(+NH+C2H2)"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c
        fragments["FA(+NH+C2H2O)"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_o
        fragments["FA(+NH+C2H4O)"] = mass - mass_o + mass_n + 5 * mass_h + 2 * mass_c + mass_o
        fragments["FA(+NH+C3H4O)"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o
        fragments["FA(+O-C2H4)"] = mass - 2 * mass_c - 4 * mass_h
        fragments["FA(+O-H2O)"] = mass - mass_h2o
        fragments["FA(+O-SL)"] = mass - mass_slfrag
        fragments["FA(+O-SL+2H)"] = mass - mass_slfrag + 2 * mass_h
        fragments["FA(+O-SL-CH2)"] = mass - mass_slfrag - mass_ch2
        fragments["FA(+O-SL-CH2-O)"] = mass - mass_slfrag - mass_ch2 - mass_o
        fragments["FA(+O-SL+O)"] = mass - mass_slfrag + mass_o
        fragments["FA(+O-SL+C2H2)"] = mass - mass_slfrag + 2 * mass_c + 2 * mass_h
        fragments["FA(+O-SL+C)"] = mass - mass_slfrag + mass_c
        fragments["FA(+O-H2O+CH2)"] = mass - mass_h2o + mass_ch2
        fragments["FA(-H)"] = acylmass - mass_h
        fragments["FA(-H2O+CH2)"] = acylmass - mass_h - mass_h2o + mass_ch2
        fragments["FA(-CH2)"] = acylmass - mass_ch2 - mass_h
        fragments["FA(-C)"] = acylmass - mass_c - mass_h
        fragments["FA(+O+NH)"] = mass + mass_n + mass_h
        fragments["FA(+O+C)"] = mass + mass_c
        fragments["FA(+2O)"] = mass + mass_o
        fragments["FA(+O+C3H6)"] = mass + 3 * mass_ch2
        fragments["FA(+O+2O-CH2)"] = mass + 3 * mass_ch2
        fragments["FA(-H+NOH2)"] = acylmass + mass_n + mass_h + mass_o
        fragments["FA(+O+O-2H)"] = mass + mass_o - 2 * mass_h

        if nhydrox > 0:
            fragments["FA(-CH2O)"] = acylmass - mass_ch2 - mass_o
            fragments["FA(+O-2H)"] = mass - 2 * mass_h

        # fragments["FA+Glycerol"] = mass + mass_glycerol - mass_h2o  # FA + glycerol - H2O
        fragments[
            "FA+GP"] = mass + mass_glycerol + mass_phosphate - 2 * mass_h2o + 2 * mass_h  # FA + glycerol + phosphate - 2 H2O + 2 H
        fragments[
            "FA+GP(-H2O)"] = mass + mass_glycerol + mass_phosphate - 3 * mass_h2o + 2 * mass_h  # FA + glycerol + phosphate - 3 H2O

        ### Neutral Losses
        nls["FA(-H)"] = mass - mass_o - mass_h  # Neutral loss of fatty acid (-H)
        nls["FA(-H)-(H2O)"] = mass + mass_h  # Neutral loss of fatty acid (-H2O)
        nls["FA+G"] = mass + mass_glycerol - mass_h2o + mass_h  # Neutral loss of FA + glycerol
        nls["FA+GP"] = fragments["FA+GP(-H2O)"] + mass_h

        # For Ceramides
        nls["FA(+O+H2O)"] = mass + mass_h2o
        nls["FA(+O+SL)"] = mass + mass_slfrag
        nls["FA(+O+NC2H3O3)"] = mass + 3 * mass_h + 2 * mass_c + mass_n + 3 * mass_o
        nls["FA(+O+CHO)"] = mass + mass_c + mass_o + mass_h

        nls["FA(+NH+C2H2+H3O)"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_h2o + mass_h
        nls["FA(+NH+C3H4O+H3O)"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o + mass_h2o + mass_h
        nls["FA(A)"] = acylmass
        nls["FA(+O)"] = mass
        nls["FA(A+OH)"] = acylmass + mass_o + mass_h
        nls["FA(A+CO2H2)"] = acylmass + 2 * mass_o + mass_c + 2 * mass_h
        nls["FA(-O)"] = mass - mass_o
        # nls["FA(-H)"] = mass - mass_o
        nls["FA(A-2CHO)"] = acylmass - mass_o - 2 * mass_c + mass_h

        if nhydrox > 0:
            nls["FA(+O+C2H3O2)"] = mass + 2 * mass_c + 2 * mass_o + 3 * mass_h
            nls["FA(+O+C2H7O2)"] = mass + 2 * mass_c + 2 * mass_o + 7 * mass_h
            nls["FA(+O+C2H5O3)"] = mass + 2 * mass_c + 3 * mass_o + 5 * mass_h
            nls["FA(+O+C2H9O3)"] = mass + 2 * mass_c + 3 * mass_o + 9 * mass_h
            nls["FA(+O+C3H7O3)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h
            nls["FA(+O+C4H10O3N1)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h
            nls["FA(+O+C2H3O2-CH2)"] = mass + 2 * mass_c + 2 * mass_o + 3 * mass_h - mass_ch2
            nls["FA(+O+C2H7O2-CH2)"] = mass + 2 * mass_c + 2 * mass_o + 7 * mass_h - mass_ch2
            nls["FA(+O+C2H5O3-CH2)"] = mass + 2 * mass_c + 3 * mass_o + 5 * mass_h - mass_ch2
            nls["FA(+O+C2H9O3-CH2)"] = mass + 2 * mass_c + 3 * mass_o + 9 * mass_h - mass_ch2
            nls["FA(+O+C3H7O3-CH2)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h - mass_ch2
            nls["FA(+O+C4H10O3N1-CH2)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h - mass_ch2
            nls["FA(+O-H)"] = mass - mass_h
            nls["FA(+O+3H)"] = mass + 3 * mass_h
            nls["FA(+O-H+H2O)"] = mass - mass_h + mass_h2o
            nls["FA(+O+3H+H2O)"] = mass + 3 * mass_h + mass_h2o
            nls["FA(+O+CH3O)"] = mass + mass_c + 3 * mass_h + mass_o
            nls["FA(+O+C2H6ON)"] = mass + 2 * mass_c + mass_n + 6 * mass_h + mass_o
            fragments[f"FA(+O-2H)"] = mass - 2 * mass_h

            # nls["FA(+O+C3H7O3-CH2)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h
            # nls["FA(+O+C4H10O3N1-CH2)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h - mass_ch2

        if adduct == "[M+CH3COO]-":
            # For Ceramides
            nls["FA(+O+H2O)+Acetate"] = mass + mass_h2o + mass_acetate + 1 * mass_h
            nls[
                "FA(+NH+C2H2+H3O)+Acetate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_h2o + mass_h + + mass_acetate + 1 * mass_h
            nls[
                "FA(+NH+C3H4O+H3O)+Acetate"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o + mass_h2o + mass_h + mass_acetate + 1 * mass_h
            nls["FA(A-2CHO)+Acetate"] = acylmass - mass_o - 2 * mass_c + mass_h + mass_acetate + mass_h

        if adduct == "[M+HCOO]-":
            nls["FA(+O+H2O)+Formate"] = mass + mass_h2o + mass_formate + 1 * mass_h
            nls[
                "FA(+NH+C2H2+H3O)+Formate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_h2o + mass_h + mass_formate + 1 * mass_h
            nls[
                "FA(+NH+C3H4O+H3O)+Formate"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o + mass_h2o + mass_h + mass_formate + 1 * mass_h
            nls["FA(A-2CHO)+Formate"] = acylmass - mass_o - 2 * mass_c + mass_h + mass_formate + mass_h

        for hgfrag in hgfrags:
            hgfrag = float(hgfrag)
            fragments[f"FA(+O)+HG({hgfrag})"] = mass + hgfrag + mass_h  # FA + headgroup fragment
            fragments[f"FA(-H)+HG({hgfrag})"] = mass - mass_o - mass_h + hgfrag  # FA + headgroup fragment

            nls[f"FA(+O)+HG({hgfrag})"] = mass + hgfrag + 2 * mass_h  # NL of FA(+O) + headgroup fragment
            nls[f"FA(-H)+HG({hgfrag})"] = mass - mass_o + hgfrag  # NL of FA(-H) + headgroup fragment

            nls[f"FA(+NH+C2H2)+HG({hgfrag})"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + hgfrag + 2 * mass_h

            if adduct == "[M+CH3COO]-":
                nls[
                    f"FA(+O)+HG({hgfrag})+Acetate"] = mass + hgfrag + mass_acetate + 2 * mass_h  # NL of FA(+O) + headgroup fragment + acetate
                nls[
                    f"FA(-H)+HG({hgfrag})+Acetate"] = mass - mass_o - mass_h + hgfrag + mass_acetate + 2 * mass_h  # NL of FA(-H) + headgroup fragment + acetate
                nls[
                    f"FA(+O)+HG({hgfrag})+Acetate+CH3"] = mass + hgfrag + mass_acetate + mass_ch2 + 2 * mass_h  # NL of FA(+O) + headgroup fragment + acetate + CH3
                nls[
                    f"FA(-H)+HG({hgfrag})+Acetate+CH3"] = mass - mass_o - mass_h + hgfrag + mass_acetate + mass_ch2 + 2 * mass_h
                nls[
                    f"FA(+NH+C2H2)+HG({hgfrag})+Acetate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + hgfrag + mass_acetate + 3 * mass_h

            if adduct == "[M+HCOO]-":
                nls[f"FA(+O)+HG({hgfrag})+Formate"] = mass + hgfrag + mass_formate + 2 * mass_h
                nls[f"FA(-H)+HG({hgfrag})+Formate"] = mass - mass_o - mass_h + hgfrag + mass_formate + 2 * mass_h
                nls[f"FA(+O)+HG({hgfrag})+Formate+CH3"] = mass + hgfrag + mass_formate + mass_ch2 + 2 * mass_h
                nls[
                    f"FA(-H)+HG({hgfrag})+Formate+CH3"] = mass - mass_o - mass_h + hgfrag + mass_formate + mass_ch2 + 2 * mass_h
                nls[
                    f"FA(+NH+C2H2)+HG({hgfrag})+Formate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + hgfrag + mass_formate + 3 * mass_h

        for hgbase in hgbases:
            hgbase = float(hgbase)
            if hgbase == 0:
                continue
            fragments[f"FA(+O)+HG({hgbase})"] = mass + hgbase  # FA + headgroup base fragment
            fragments[f"FA(-H)+HG({hgbase})"] = mass - mass_h2o + hgbase  # FA + headgroup base fragment

            nls[f"FA(+O)+HG({hgbase})"] = mass + hgbase + mass_h  # NL of FA(+O) + headgroup base fragment
            nls[f"FA(-H)+HG({hgbase})"] = mass - mass_h2o + hgbase + mass_h  # NL of FA(-H) + headgroup base fragment

            if adduct == "[M+CH3COO]-":
                nls[
                    f"FA(+O)+HG({hgbase})+Acetate"] = mass + hgbase + mass_acetate + 2 * mass_h  # NL of FA(+O) + headgroup fragment + acetate
                nls[
                    f"FA(-H)+HG({hgbase})+Acetate"] = mass - mass_o - mass_h + hgbase + mass_acetate + 2 * mass_h  # NL of FA(-H) + headgroup fragment + acetate
                nls[
                    f"FA(+O)+HG({hgbase})+Acetate+CH3"] = mass + hgbase + mass_acetate + mass_ch2 + 2 * mass_h  # NL of FA(+O) + headgroup fragment + acetate + CH3
                nls[
                    f"FA(-H)+HG({hgbase})+Acetate+CH3"] = mass - mass_o - mass_h + hgbase + mass_acetate + mass_ch2 + 2 * mass_h
            if adduct == "[M+HCOO]-":
                nls[
                    f"FA(+O)+HG({hgbase})+Formate"] = mass + hgbase + mass_formate + 2 * mass_h  # NL of FA(+O) + headgroup fragment + formate
                nls[
                    f"FA(-H)+HG({hgbase})+Formate"] = mass - mass_o - mass_h + hgbase + mass_formate + 2 * mass_h  # NL of FA(-H) + headgroup fragment + formate
                nls[
                    f"FA(+O)+HG({hgbase})+Formate+CH3"] = mass + hgbase + mass_formate + mass_ch2 + 2 * mass_h  # NL of FA(+O) + headgroup fragment + formate + CH3
                nls[
                    f"FA(-H)+HG({hgbase})+Formate+CH3"] = mass - mass_o - mass_h + hgbase + mass_formate + mass_ch2 + 2 * mass_h  # NL of FA(-H) + headgroup fragment + formate + CH3

        if adduct == "[M+CH3COO]-":
            nls["FA(-H)+Acetate"] = mass - mass_h2o + mass_acetate + 2 * mass_h  # Neutral loss of FA + acetate
            nls["FA(+O)+Acetate"] = mass + mass_acetate + 2 * mass_h  # Neutral loss of FA(+O) + acetate
            nls[
                "FA(-H)+Acetate+CH3"] = mass - mass_h2o + mass_acetate + mass_ch2 + 2 * mass_h  # Neutral loss of FA + acetate + CH3
            nls[
                "FA(+O)+Acetate+CH3"] = mass + mass_acetate + mass_ch2 + 2 * mass_h  # Neutral loss of FA(+O) + acetate + CH3

        if adduct == "[M+HCOO]-":
            nls["FA(-H)+Formate"] = mass - mass_h2o + mass_formate + 2 * mass_h  # Neutral loss of FA + formate
            nls["FA(+O)+Formate"] = mass + mass_formate + 2 * mass_h  # Neutral loss of FA(+O) + formate

            nls[
                "FA(-H)+Formate+CH3"] = mass - mass_h2o + mass_formate + mass_ch2 + 2 * mass_h  # Neutral loss of FA + formate + CH3
            nls[
                "FA(+O)+Formate+CH3"] = mass + mass_formate + mass_ch2 + 2 * mass_h  # Neutral loss of FA(+O) + formate + CH3

    # Round all computed fragment and neutral loss values once here
    for key in list(fragments.keys()):
        fragments[key] = round(fragments[key], 4)
    for key in list(nls.keys()):
        nls[key] = round(nls[key], 4)

    return fragments, nls

def gen_fragdict(cname, adduct, tails, hgdf):
    # print(cname, adduct)
    if cname in hgbasedict:
        hgbase = hgbasedict[cname]
    else:
        hgbase = 0

    # Unpack headgroup fragments from hgdf
    hgrow = hgdf[(hgdf["CompoundClass"] == cname) & (hgdf["Adduct"] == adduct)]
    if len(hgrow) == 0:
        # print("No headgroup fragments found for", cname, adduct)
        hgfrags = []
        hgnls = []
    else:
        hgfrags = hgrow.iloc[0]["HeadgroupFragment_mz"]
        hgnls = hgrow.iloc[0]["HeadgroupNL_mz"]
        # Drop nans
        if pd.isna(hgfrags) or hgfrags.strip() == "":
            hgfrags = []
        else:
            hgfrags = [float(f) for f in hgfrags.split(",")]
        if pd.isna(hgnls) or hgnls.strip() == "":
            hgnls = []
        else:
            hgnls = [float(nl) for nl in hgnls.split(",")]

    # Convert to dict with "HG(FRAG1)": mz format
    fragdict = {}
    nldict = {}
    for i, frag in enumerate(hgfrags):
        fragdict[f"HG({frag})"] = frag
    for i, nl in enumerate(hgnls):
        nldict[f"HG_NL({nl})"] = nl

    # Merge these into the main dicts with tail info
    tname = ["T1", "T2", "T3", "T4"]
    tailfrags = []
    tailnls = []
    for i, tail in enumerate(tails):
        if pd.isna(tail) or tail.strip() == "":
            tailfrags.append([])
            tailnls.append([])
            continue
        predicted_frags, predicted_nls = predict_tail_fragments(tail, hgfrags, hgbase, adduct=adduct, classname=cname)
        tailfrags.append(predicted_frags)
        tailnls.append(predicted_nls)

        # Add headgroup frags/nls to fragdict with appended tail name
        for frag_name, frag_mz in predicted_frags.items():
            fragdict[f"{tname[i]}_{tail}_{frag_name}"] = frag_mz
        for nl_name, nlval in predicted_nls.items():
            nldict[f"{tname[i]}_{tail}_NL[{nl_name}]"] = nlval
    return fragdict, nldict


def assign_fragments(row, tol=0.05):
    tail = row["Tail"]
    fraglist = row["TailFragment_mz"]
    nllist = row["TailNL_mz"]
    hgfrags = row["HeadgroupFragment_mz"]
    cname = row["CompoundClass"]
    adduct = row["Adduct"]
    ignorefrag = ignore_fragments[cname] if cname in ignore_fragments else []
    ignorenl = ignore_nls[cname] if cname in ignore_nls else []
    # print(cname, adduct)
    if cname in hgbasedict:
        hgbase = hgbasedict[cname]
    else:
        hgbase = 0

    predicted_frags, predicted_nls = predict_tail_fragments(tail, hgfrags, hgbase, adduct=adduct, classname=cname)

    if pd.isna(fraglist) or fraglist.strip() == "":
        outstring = ""
    else:
        fraglist = [float(f) for f in fraglist.split(",")]

        outfraglist = []
        for frag in fraglist:
            if frag in ignorefrag:
                continue

            found = False
            for frag_name, frag_mz in predicted_frags.items():
                if abs(frag - frag_mz) <= tol:
                    outfraglist.append(frag_name)
                    found = True
            if not found:
                outfraglist.append("NOTFOUND(" + str(frag) + "|" + tail + ")")
                print("NOTFOUND Fragment:", frag, "for", cname, "Tail:", tail)
                defect = frag - predicted_frags['FA(+O)']
                print(defect)
                # form = calc_fromula_from_mass([defect], searchspace=reallysmallsearch, tolerance=200)
                # print(form)
                print(frag - predicted_frags['FA(+O)'], ",")

        outstring = ",".join(outfraglist)

    if pd.isna(nllist) or nllist.strip() == "":
        outstring2 = ""
    else:
        nllist = [float(nl) for nl in nllist.split(",")]

        outnllist = []
        for nl in nllist:
            if nl in ignorenl:
                # print("Ignoring NL:", nl, "for", cname)
                continue
            found = False
            for nl_name, nl_mz in predicted_nls.items():
                if abs(nl - nl_mz) <= tol:
                    outnllist.append(nl_name)
                    found = True
            if not found:
                outnllist.append("NOTFOUND(" + str(nl) + "|" + tail + ")")
                print("NOTFOUND Neutral Loss:", nl, "for", cname, "Tail:", tail)
                defect = nl - predicted_nls['FA(+O)']
                print(defect)
                # form = calc_fromula_from_mass([defect], searchspace=reallysmallsearch, tolerance=200)
                # print(form)
                print(nl - predicted_frags['FA(+O)'] - mass_h, ",")

        outstring2 = ",".join(outnllist)

    row["TailFragment_Assigned"] = outstring
    row["TailNL_Assigned"] = outstring2
    return row


def assign_fragments_to_df(df, tol=0.05):
    print("Assigning:", len(df))
    # print("Unique Classes:", df["CompoundClass"].unique())
    # df = df.apply(assign_fragments, axis=1, tol=tol)
    rows = []
    for index, row in df.iterrows():
        newrow = assign_fragments(row, tol=tol)
        rows.append(newrow)
    df = pd.DataFrame(rows)
    return df


def merge_fragments_df(df):
    print("Merging:", len(df))
    # Merge all with the same CompoundClass, Adduct, IonMode, HeadgroupFragment_mz,	HeadgroupNL_mz,	TailFragment_Assigned,	TailNL_Assigned
    # Drop columns TailFragment_mz, TailNL_mz, Tail
    df = df.drop(columns=["TailFragment_mz", "TailNL_mz", "Tail"])
    # print("Unique Classes:", df["CompoundClass"].unique())
    df = df.drop_duplicates()
    # print("Unique Classes:", df["CompoundClass"].unique())
    print("Finished:", len(df))
    return df


if __name__ == "__main__":
    import time

    starttime = time.perf_counter()
    os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

    df = pd.read_csv("LipidTailFragments.csv", dtype=str)

    # Select df for negative ion mode
    df = df[df["IonMode"] == "negative"].reset_index(drop=True)
    # Select df for specific class
    # df = df[df["CompoundClass"] == "SL"].reset_index(drop=True)
    # Isolate tails with 18:1
    # df = df[df["Tail"] == "18:1;O"].reset_index(drop=True)

    df = assign_fragments_to_df(df, tol=0.05)
    df = merge_fragments_df(df)
    # Drop rows where both TailFragment_Assigned and TailNL_Assigned are empty
    # df = df[
    #     ~((df["TailFragment_Assigned"].isna() | (df["TailFragment_Assigned"].str.strip() == "")))].reset_index(
    #     drop=True)
    # print(df)
    # Write to xlsx
    # df.to_excel("LipidTailFragments_18120.xlsx", index=False)
    # print("Written output file.")
    # Send df to list and print

    writer = sf.ExcelWriter("LipidTailFragments_Test.xlsx")
    sfdf = sf.StyleFrame(df)
    # df.to_excel(writer, sheet_name='Sheet 1', index=False, na_rep='NaN')
    sfdf.to_excel(writer, best_fit=df.columns.tolist())

    writer.close()

    print(time.perf_counter() - starttime)
