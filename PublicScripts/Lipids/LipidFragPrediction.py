import pandas as pd
import numpy as np


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

hgbasedict = {
    "Ac2PIM1": [324.10226, 460.1001],
    "AHexCer": [208.052, 222.08, 2.96706, 145.0506, 163.06144, 181.07144],
    "ASM": [4.03, 4.04, 6.04, 6.05, 8.06, 8.07, 10.08, 12.09, 12.1, 62.03, 62.04, 64.05,
            64.06, 66.06, 66.07, 68.08, 68.09, 70.1, 72.11, 72.12, 76.05, 76.06, 78.06, 78.07,
            80.08, 80.09, 82.1, 84.11, 84.12, 86.13, 86.14, 162.0261],
    "BMP": [172.04144],
    "CerP": [98.98144],
    "Cer_AP": [89.04366, 135.05366, 71.0315, -53.0281, 1.0119, -1.01016, 16.999834, 35.00984],
    "Cer_AS": [28.97706, 103.02366, 117.04366, 89.04366, 57.02366, 42.98975],
    "Cer_BDS": [9.9848, 34.00206, 11.99994, -44.02434, -12.00434],
    "Cer_BS": [93.04246, 46.00116, -14.01724, 33.02276, -44.03154],
    "Cer_EBDS": [73.0528],
    "Cer_EOS": [39.997, 33.98204, 45.98204],
    "Cer_EODS": [39.997],
    "Cer_NDS": [31.02144],
    "Cer_NS": [31.02144],
    "CL": [71.0111, 340.01894],
    "DGCC": [234.1336],
    "DGDG": [342.14144],
    "DGGA": [250.07, 194.06144],
    "DGTS": [218.1336],
    "DLCL": [210.03, 364.032],
    "EtherDGDG": [426.13876, 380.1323, 398.1423, 342.14144, 163.06144, 241.0726, 403.1226],
    "EtherMGDG": [236.08660, 264.083060, 278.10016, 180.08794, 241.0661],
    "EtherPE": [197.04, 179.03016],
    "EtherPG": [228.04 - mass_o, 194.03206, 210.02856],
    "EtherPI": [316.05 - mass_o, 316.05, 162.05],
    "EtherPS": [223.02716, 87.032 - mass_o, 223.0022 - mass_o, 87.032],
    "EtherSMGDG": [298.04069],
    "EtherOxPC": [-36.0225],
    "EtherOxPE": [-36.0225, 197.0464, 179.03856],
    "GM3": [206.94144, 651.24144],
    "Hex2Cer": [338.0849, 176.0349, 294.05706, 204.02706, 308.07706, 218.04706, 343.12144, 355.12144],
    "Hex3Cer": [500.1349, 338.0849, 176.0349, 366.07706, 380.09706, 505.17143, 517.17144],
    "HexCer": [176.0349],
    "HexCer_AP": [71.0315, 212.05066, 182.04066, 237.08066, 283.09066, 149.02216, 179.03216, 204.06216,
                  250.07216, 190.05216, 144.04216, 119.01216, 89.00216, 236.05216, 215.07443, -53.0206],
    "HexCer_EOS": [193.04706, 207.06706, 12.9895, 191.0443, 158.0197, 176.0298, 204.02196, 190.05066, 208.03804,
                   178.02374, 196.03064],
    "HexCer_HDS": [193.04706, 207.06706, 12.9895, 196.02974, 208.02974],
    "HexCer_HS": [193.04706, 207.06706, 12.9895, 196.02974, 208.02974],
    "HexCer_NDS": [193.07144, 181.07144, 163.06144],
    "HexCer_NS": [193.07144, 181.07144, 163.06144, 146.03264],
    "LNAPE": [179.0326, 197.0483],
    "LNAPS": [69.02086],
    "MGDG": [218.077, 180.09144],
    "NAGly": [43.99206, 94.04854, 111.07444],
    "NAGlySer": [61.99796, 30.00796, 151.06853, 180.09854, 36.04414, 123.07573, 123.07734, 123.07874, 123.07444,
                 123.08014],
    "NAOrn": [133.10014],
    "OxPC": [-36.0225],
    "OxPE": [-36.0225],
    "OxPI": [-36.0225],
    "OxPS": [87.034, 223.0022, -36.0225],
    "OxPG": [228.04, -36.0225],
    "PC": [281.102, 240.1036, 71.0057, 92.9958, 216.9858, 130.08363, 190.97364],
    "PE": [198.05144],
    "PEtOH": [165.0336, 200.07144, 126.03144],
    "PE_Cer": [-3.04804],
    "PI": [316.05],
    "PI_Cer": [-3.04804],
    "PG": [228.04, 172.04144],
    "PMeOH": [151.0136, 126.03144, 186.05144],
    "PS": [87.032, 223.0022, 242.0436],
    "SHexCer": [273.03144, 261.03144, 243.02144, 80.96144, 275.98974, 293.01974, 225.98674, 145.0445],
    "SL": [79.9594, 110.0236, 62.9479],
    "SM": [184.07144],
    "SQDG": [282.04, 283.04356, 318.08144, 244.05144],
}

ignore_nls = {"Hex2Cer": [46.01, 46.0], "Hex3Cer": [46.01], "HexCer": [46.01], "AHexCer": [46.01], "SM": [60.02],
              "SL": [111.98], "DG": [17.03, 35.04], "SHexCer": [79.96, 242.01], "TG": [17.03], "GM3": [488.19]}
acelated_classes = ["AHexCer", "ASM"]
adduct_nls = {"[M+CH3COO]-": [60.02], "[M+HCOO]-": [46.01], "[M+NH4]+": [17.0266, 35.037]}


def predict_tail_fragments(tail, hgfrags=[], hgbases=[], hgnls=[], adduct="[M-H]-", classname="", mode="negative"):
    # If hgfrags is a string, split by comma
    if isinstance(hgfrags, str):
        hgfrags = hgfrags.split(",")
    if not isinstance(hgfrags, list):
        hgfrags = [hgfrags]

    if isinstance(hgnls, str):
        hgnls = hgnls.split(",")
    if not isinstance(hgnls, list):
        hgnls = [hgnls]

    if classname != "":
        hgbases = hgbasedict[classname] if classname in hgbasedict else []
    else:
        hgbases = []

    if isinstance(hgbases, str):
        hgbases = hgbases.split(",")
    if not isinstance(hgbases, list):
        hgbases = [hgbases]

    if hgbases == [] and classname in hgbasedict:
        hgbases = hgbasedict[classname]

    fragments = {}
    nls = {}
    if pd.isna(tail):
        return fragments, nls
    if ":" not in tail:
        # print("Error: No ':' in tail string:", tail)
        return fragments, nls

    parts = tail.split(":")
    numox = 0
    nhydrox = 0
    nmethyl = 0
    ndeut = 0
    nbound = False
    obound = False
    pbound = False
    sn = ""
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

    if "P-" in parts[0]:
        parts[0] = parts[0].replace("P-", "")
        pbound = True

    if "N-" in parts[0]:
        parts[0] = parts[0].replace("N-", "")
        nbound = True

    if "-SN1" in parts[1]:
        parts[1] = parts[1].replace("-SN1", "")
        sn = "SN1"
    elif "-SN2" in parts[1]:
        parts[1] = parts[1].replace("-SN2", "")
        sn = "SN2"

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
    if nmethyl > 0:
        mass += nmethyl * (mass_c + 2 * mass_h)  # Add methyl group

    # If O bound, remove O and add 2 H
    if obound and classname not in acelated_classes:
        mass += - mass_o + 2 * mass_h

    if pbound:
        mass += -mass_o

    # print(numox, nhydrox, nmethyl, ndeut, num_carbons, num_double_bonds)

    if mode == "negative":
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


    elif mode == "positive":
        fragments["FA(+OH)"] = mass + mass_h
        fragments["FA(-H)"] = mass - mass_h2o + mass_h
        fragments["FA(+OH-OH)"] = mass - mass_h2o + 2 * mass_h
        fragments["FA(+G)"] = mass + mass_glycerol - mass_h2o + mass_h
        fragments["FA(+G-OH)"] = mass + mass_glycerol - 2 * mass_h2o + 2 * mass_h
        fragments["FA(+GP)"] = mass + mass_glycerol + mass_phosphate - 2 * mass_h2o + 3 * mass_h
        fragments["FA(+GP-H2O)"] = mass + mass_glycerol + mass_phosphate - 3 * mass_h2o + 3 * mass_h
        fragments["FA(+O-O2H2)"] = mass - 2 * mass_o - 2 * mass_h
        fragments["FA(+O-CO2H2)"] = mass - 2 * mass_o - 2 * mass_h - mass_c
        fragments["FA(+O-C2NOH5)"] = mass - 2 * mass_c - mass_n - 5 * mass_h - mass_o
        fragments["FA(+O+C2H8NO3P)"] = mass + 2 * mass_c + mass_n + 8 * mass_h + mass_p + 3 * mass_o
        fragments["FA(+H+C2H4N)"] = mass + 2 * mass_ch2 + mass_n - mass_o + mass_h

        nls["FA(+O)"] = mass
        nls["FA(+OH)"] = mass + mass_h
        nls["FA(+H)"] = mass - mass_o
        nls["FA(-H)"] = mass - mass_h2o + mass_h
        nls["FA(+O+H2O)"] = mass + mass_h2o
        nls["FA(+OH+H2O)"] = mass + mass_h2o + mass_h

        for hgbase in hgbases:
            hgbase = float(hgbase)
            if hgbase == 0:
                continue
            fragments[f"FA(+OH)+HG({hgbase})"] = mass + hgbase + mass_h
            fragments[f"FA(-H)+HG({hgbase})"] = mass - mass_h2o + hgbase + mass_h

            nls[f"FA(+OH)+HG({hgbase})"] = mass + hgbase
            nls[f"FA(-H)+HG({hgbase})"] = mass - mass_h2o + hgbase
            nls[f"FA(+OH-CH2)+HG({hgbase})"] = mass + hgbase - mass_ch2

    else:
        print("Error: mode must be 'negative' or 'positive'")
        return fragments, nls

    # Round all computed fragment and neutral loss values once here
    for key in list(fragments.keys()):
        fragments[key] = round(fragments[key], 4)
    for key in list(nls.keys()):
        nls[key] = round(nls[key], 4)

    return fragments, nls


def gen_fragdict(cname, adduct, tails, hgdf, mode=None):
    if mode is None:
        if "+" in adduct:
            mode = "positive"
        else:
            mode = "negative"
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
        predicted_frags, predicted_nls = predict_tail_fragments(tail, hgfrags, hgbase, adduct=adduct, classname=cname, mode=mode)
        tailfrags.append(predicted_frags)
        tailnls.append(predicted_nls)

        # Add headgroup frags/nls to fragdict with appended tail name
        for frag_name, frag_mz in predicted_frags.items():
            fragdict[f"{tname[i]}_{tail}_{frag_name}"] = frag_mz
        for nl_name, nlval in predicted_nls.items():
            nldict[f"{tname[i]}_{tail}_NL[{nl_name}]"] = nlval
    return fragdict, nldict

def merge_lists(fraglist, nllist):
    # Remove NOTFOUND entries where they were found in the other list
    if len(fraglist) != len(nllist):
        print("Error: fraglist and nllist have different lengths")
        print(len(fraglist), len(nllist))
        raise ValueError

    for i in range(len(fraglist)):
        frag = fraglist[i]
        nl = nllist[i]
        if "NOTFOUND" in frag and not ("NOTFOUND" in nl):
            fraglist[i] = nl
        elif "NOTFOUND" in nl and not ("NOTFOUND" in frag):
            nllist[i] = frag
    return fraglist, nllist



def assign_df(df, hgdf):
    if "T1" not in df.columns:
        df = set_basic_tail_names(df, columnname="Metabolite name")

    outrows = []
    for i, row in df.iterrows():
        row = assign_ref_match(row, tol=0.05, hgdf=hgdf)
        outrows.append(row)

    return pd.DataFrame(outrows)

def assign_ref_match(row, tol=0.05, hgdf=None, verbose=True):
    # Unpack row
    name = row["Metabolite name"]
    cname = row["Ontology"]
    adduct = row["Adduct type"]
    tail1 = row["T1"]
    tail2 = row["T2"]
    tail3 = row["T3"]
    tail4 = row["T4"]
    tails = [tail1, tail2, tail3, tail4]
    msms_mz = row["MS/MS spectrum"]
    ref_spectrum = row["Ref Spec"]
    ref_mz = float(row["Reference m/z"])

    # Unpack reference spectrum and msms_mz
    mzlist = [m.split(":") for m in msms_mz.split(" ") if m.strip() != ""]
    reflist = [m.split(":") for m in ref_spectrum.split(" ") if m.strip() != ""]
    mzlist = np.array([[float(mz[0]), float(mz[1])] for mz in mzlist])
    reflist = np.array([[float(r[0]), float(r[1])] for r in reflist])
    # Filter each to 1% intensity threshold
    mzthreshold = np.max(mzlist[:, 1]) * 0.01
    refthreshold = np.max(reflist[:, 1]) * 0.01
    mzlist = mzlist[mzlist[:, 1] >= mzthreshold]
    reflist = reflist[reflist[:, 1] >= refthreshold]

    matchedmz = []
    for r in reflist[:, 0]:
        for m in mzlist[:, 0]:
            if check_match(m, r, tol):
                matchedmz.append(r)
                break
    matchedmz = np.array(matchedmz)
    # Filter to >5 m/z below precursor
    matchedmz = matchedmz[matchedmz < ref_mz - 5]
    matchednl = ref_mz - matchedmz

    # Get the potential adduct neutral losses
    adductnls = adduct_nls[adduct] if adduct in adduct_nls else []

    fragdict, nldict = gen_fragdict(cname, adduct, tails, hgdf)

    # if verbose:
    #     print(name, fragdict, nldict)

    nfrag_missed = 0
    nnl_missed = 0
    tailcounts = np.zeros(4, dtype=int)
    hgcounts = 0
    adductcounts = 0

    outfraglist = []
    for frag in matchedmz:
        found = False
        for frag_name, frag_mz in fragdict.items():
            if found:
                continue
            if check_match(frag, frag_mz, tol):
                outfraglist.append(frag_name)
                found = True
                if frag_name[:2] == "T1":
                    tailcounts[0] += 1
                elif frag_name[:2] == "T2":
                    tailcounts[1] += 1
                elif frag_name[:2] == "T3":
                    tailcounts[2] += 1
                elif frag_name[:2] == "T4":
                    tailcounts[3] += 1
                elif frag_name[:2] == "HG":
                    hgcounts += 1
                else:
                    print("Unknown fragment name:", frag_name)
        if not found:
            outfraglist.append("NOTFOUND[" + str(frag) + "|" + name + "]")
            nfrag_missed += 1

    tailnlcounts = np.zeros(4, dtype=int)
    hgcountsnl = 0
    adductcountsnl = 0
    outnllist = []

    for nl in matchednl:
        found = False
        for nl_mz in adductnls:
            if found:
                continue
            if check_match(nl, nl_mz, tol):
                outnllist.append("Adduct[" + str(nl) + "]")
                adductcountsnl += 1
                found = True
                continue

        for nl_name, nl_mz in nldict.items():
            if found:
                continue
            if check_match(nl, nl_mz, tol):
                outnllist.append(nl_name)
                found = True
                if nl_name[:2] == "T1":
                    tailnlcounts[0] += 1
                elif nl_name[:2] == "T2":
                    tailnlcounts[1] += 1
                elif nl_name[:2] == "T3":
                    tailnlcounts[2] += 1
                elif nl_name[:2] == "T4":
                    tailnlcounts[3] += 1
                elif nl_name[:2] == "HG":
                    hgcountsnl += 1
                else:
                    print("Unknown NL name:", nl_name)
        if not found:
            outnllist.append("NOTFOUND[" + str(nl) + "|" + name + "]")
            nnl_missed += 1

    outfraglist, outnllist = merge_lists(outfraglist, outnllist)

    if verbose:
        for frag in outfraglist:
            if "NOTFOUND" in frag:
                print("NOTFOUND Fragment:", frag, "for", cname, "Name:", name, "Adduct:", adduct)
        for nl in outnllist:
            if "NOTFOUND" in nl:
                print("NOTFOUND Neutral Loss:", nl, "for", cname, "Name:", name, "Adduct:", adduct)

    outstring = ",".join(outfraglist)
    outstring2 = ",".join(outnllist)

    for i in range(4):
        row[f"T{i+1}_Fragments_Found"] = tailcounts[i]
        row[f"T{i+1}_NL_Found"] = tailnlcounts[i]
    row["HG_Fragments_Found"] = hgcounts
    row["HG_NL_Found"] = hgcountsnl
    row["Adduct_NL_Found"] = adductcountsnl
    row["Fragments_Assigned"] = outstring
    row["NL_Assigned"] = outstring2
    row["Fragments_Missed"] = nfrag_missed
    row["NL_Missed"] = nnl_missed

    hgmatch = 0
    if hgcounts > 0 or hgcountsnl > 0:
        hgmatch = 1
    row["Headgroup_Match"] = hgmatch

    tailmatch = 0
    ntails = 0
    for idx, t in enumerate(tails):
        if t.strip() == "":
            continue
        ntails += 1
        if tailcounts[idx] > 0 or tailnlcounts[idx] > 0:
            tailmatch += 1
    row["Tail_Match_Count"] = tailmatch
    row["Total_Tails"] = ntails
    row["Tail_Match_Fraction"] = tailmatch / ntails if ntails > 0 else 0

    adduct_match = 0
    if adductcountsnl > 0:
        adduct_match = 1
    row["Adduct_Match"] = adduct_match

    grade = "F: No Match"
    if hgmatch == 1 and tailmatch == ntails and adduct_match == 1:
        grade = "A+: Full Match"
    elif hgmatch == 1 and tailmatch == ntails:
        grade = "A: Headgroup and All Tails Match"
    elif hgmatch == 1 and tailmatch > 0:
        grade = "B: Headgroup and Some Tails Match"
    elif hgmatch == 1:
        grade = "C: Headgroup Match Only"
    elif tailmatch == ntails and ntails > 0:
        grade = "C: All Tails Match Only"
    elif tailmatch > 0 and ntails > 0:
        grade = "D: Some Tails Match Only"
    row["Match_Grade"] = grade
    return row

def check_match(query_mz, target, tol=0.05):
    if abs(query_mz - target) <= tol:
        return True
    return False