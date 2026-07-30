import numpy as np
import os
import re
import pandas as pd
import lxml.etree
import matplotlib.pyplot as plt
from sklearn.covariance import MinCovDet
from scipy.spatial.distance import mahalanobis
import networkx as nx
import matchms as mms
import matplotlib as mpl
import mplcursors
from unidec.tools import center_of_mass
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import rdkit.Chem.MACCSkeys as MACCSkeys
import time

# import warnings
# warnings.filterwarnings("error")



class_color_map = {
    "PC": "#788BFF",
    "OxPC": "#788BFF",
    "PC O": "#BFD7FF",
    "EtherPC": "#BFD7FF",
    "LPC-O/P": "#27D3F5",
    "LPC": "#3019E3",
    "EtherLPC": "#27D3F5",
    "PE": "#6017A9",
    "PE-O": "#A276E5",
    "EtherPE": "#A276E5",
    "LPE": "#8F5CD2",
    "EtherLPE": "#8F5CD2",
    "PEtOH" : "#D291EF",
    "PMeOH": "#D291EF",
    "PG": "#356935",
    "LPG": "#B1D8B1",
    "LPS": "#FDF32F",
    "PS": "#C89504",
    "PA": "#95F4D4",
    "LPA": "#42EBB3",
    "PI": "#8FFBFF",
    "EtherPI": "#8FFBFF",
    "EtherDG": "#F4845F",
    "EtherPS": "#C89504",
    "LPI": "#00CAD1",
    "CL": "#FFA200",
    "MLCL": "#FFDAA5",
    "DLCL": "#FFE6C2",
    "Cer": "#AB30A1",
    "Cer_NS": "#AB30A1",
    "Cer_NDS": "#AB30A1",
    "Cer_EODS": "#AB30A1",
    "Cer_EOS": "#AB30A1",
    "Cer_EBDS": "#AB30A1",
    "Cer_AP": "#AB30A1",
    "Cer_AS": "#AB30A1",
    "CAR": "#AB30A1",
    "HexCer": "#F708E3",
    "HexCer_NS": "#F708E3",
    "HexCer_NDS": "#F708E3",
    "HexCer_HDS": "#F708E3",
    "HexCer_AP": "#F708E3",
    "HexCer_HS": "#F708E3",
    "Hex2Cer": "#FF73FF",
    "SM": "#FFDB73",
    "ASM": "#FFDB73",
    "CoQ": "#F25C54",
    "MAG": "#F7B267",
    "DG": "#F4845F",
    "DAG": "#F4845F",
    "TG": "#8F1458",
    "TAG": "#8F1458",
    "FA": "#FFA1DB",
    "FAHFA": "#FF0000",
    "CE": "#D10081",
    "NAE": "#777777",
    "ST": "#B40808",
    "Standard": "#0D0D0D",
    "Other": "#333333",
    "GM3": "#FF73FF",
    "HBMP": "#FFA200",
    "SHexCer": "#F708E3",
    "OxPI": "#8FFBFF",
    "PI_Cer": "#8FFBFF",
    "PE_Cer": "#6017A9",
    "AHexBRS": "#FF73FF",
    "AHexCAS": "#FF73FF",
    "EtherPG": "#356935",
    "OxPE": "#6017A9",
    "AHexCer": "#FF73FF",
    "EtherLPG": "#B1D8B1",
    "EtherPG": "#356935",
    "EtherTG": "#8F1458",
    "BileAcid": "#FF8800",
    # Red
    "SL": "#B40808",
    "LNAPE": "#8F5CD2",
    "LNAPS": "#FDF32F",
    "EtherDGDG": "#FFA200",
    "EtherMGDG": "#FFA200",
    "DGGA": "#FFA200",
    "SQDG": "#FFA200",
    "MGDG": "#FFA200",
    "OxPS": "#C89504",
    "DGDG" : "#FFA200",
    "SHex": "#F708E3",
    "Sph": "#AB30A1",
    "ADGGA": "#FFA200",
    "STSE" : "#B40808",
    "GL": "#F7B267",
}

# python
default_color = "#0D0D0D"
lower_color_map = {k.strip().lower(): v for k, v in class_color_map.items()}

def get_color_for_class(cls):
    key = str(cls).strip().lower()
    color = lower_color_map.get(key)
    if color is None:
        print(f"WARNING: No color found for class '{cls}'. Using default color.")
        return default_color
    return color

posadducts = [
    ["[M+H]+", 1, 1.00782503207, True],
    ["[M+NH4]+", 1, 18.03437413, True],
    ["[M+Na]+", 1, 22.9897692809, True],
    ["[M+CH3OH+H]+", 1, 33.03403978207, True],
    ["[M+K]+", 1, 38.9637066, True],
    ["[M+Li]+", 1, 7.01604545, False],
    ["[M+ACN+H]+", 1, 42.03473142507, False],
    ["[M+H2O]+", 1, 17.00273964793, False],
    ["[M+Zn+H]+", 1, 65.919582379, False],
    ["[M+IsoProp+Na+H]+", 1, 81.0596798687, False],
    ["[M+C6H100S+H]+", 1, -161.0449838463, False],
    ["[2M+H]+", 1, 1.00782503207, True],
    ["[2M+NH4]+", 1, 18.03437413, True],
    ["[2M+Na]+", 1, 22.9897692809, True],
    ["[2M+ACN+H]+", 1, 42.0347312307, False],
    ["[2M+ACN+Na]+", 1, 64.061638839, False],
    ["[M+2H]+", 2, 2.01565006414, False],
    ["[M+H+Na]+", 2, 19.04219961207, False],
    ["[M+H+K]+", 2, 39.971517127, False],
    ["[M+2ACN+H]+", 2, 43.04219191614, False],
    ["[M+2ACN+Na]+", 2, 64.0674826414, False],
    ["[M+3ACN+H]+", 3, 125.02597294614, False],
    ["[M+3H]+", 3, 3.02347509261, False],
    ["[M+3Na]+", 3, 68.9693078427, False],
]

negadducts = [
    ["[M-H]-", 1, -1.00782503207, True],
    ["[M-H2O-H]-", 1, -19.01838971207, False],
    ["[M+Na-2H]-", 1, 20.97411921676, False],
    ["[M+Cl]-", 1, 34.96885268, False],
    ["[M+K-2H]-", 1, 36.94805661586, False],
    ["[M+HCOO]-", 1, 44.997654, False],
    ["[M+CH3COO]-", 1, 59.013305, True],
    ["[M+C2H3N+Na-2H]-", 1, 62.00066831777, True],
    ["[M+Br]-", 1, 78.9183371, False],
    ["[M+TFA]-", 1, 112.98503896793, False],
    ["[M-C6H10O4-H]-", 1, -147.06573383101, False],
    ["[M-C6H10O5-H]-", 1, -163.06064845057, False],
    ["[M-C6H8O6-H]-", 1, -177.03991305059, False],
    ["[M+CH3COONa-H]-", 1, 80.99524996793, False],
    ["[2M-H]-", 1, -1.00782503207, True],
    ["[2M+FA-H]-", 1, 44.99765496793, False],
    ["[2M+Hac-H]-", 1, 59.01330596793, False],
    ["[3M-H]-", 1, -1.00782503207, True],
    ["[M-2H]2-", 2, -2.01565006414, False],
    ["[M-3H]3-", 3, -3.02347509261, False],
]

adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M+HCOO]-", "[M+HCOO]", "-1", -1, "[M+HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"],
                              ["[M+CH3COO]-", "[M+CH3COO]", "-1", -1, "[M+CH3COO]1-"],
                              ["[M+Na]+", "[M+Na]", "1", -1, "[M+Na]1+"],
                              ["[M+H-H2O]+", "[M+H-H2O]", "+1", 1, "[M+H-H2O]1+"],
                              ["[M7H2+H-H2O]+", "[M7H2+H-H2O]", "+1", 1, "[M7H2+H-H2O]1+"]
                              ])

class_nfattyacids = {
    "CL": 4,
    "Cer": 2,
    "Ch": 0,
    "ChE": 1,
    "DAG": 2,
    "HexCer": 2,
    "LPA": 1,
    "LPC": 1,
    "LPC-O/P": 1,
    "LPE": 1,
    "LPE-O/P": 1,
    "LPG": 1,
    "LPI": 1,
    "LPS": 1,
    "LSM": 1,
    "MAG": 1,
    "PA": 2,
    "PC": 2,
    "PC-O/P": 2,
    "PE": 2,
    "PE-O/P": 2,
    "PG": 2,
    "PI": 2,
    "PS": 2,
    "SHexCer": 2,
    "SM": 2,
    "ST": 1,
    "TAG": 3
}

hg_simplifier = {"Cer_AP": "Cer", "Cer_HDS": "Cer", 'Cer_NDS': "Cer", 'Cer_NP': "Cer", 'Cer_NS': "Cer",
                 "Cer_EODS" : "Cer", "Cer_EOS" : "Cer", "Cer_HS" : "Cer", "Cer_AS": "Cer", "Cer_ADS" : "Cer",
                 "HexCer_AP": "HexCer", "HexCer_NDS": "HexCer", "HexCer_HS": "HexCer", "HexCer_NS": "HexCer",
                 "HexCer_EOS": "HexCer", "HexCer_HDS": "HexCer", "AHexCer": "HexCer", "Hex2Cer": "HexCer",
                 "Hex3Cer": "HexCer", "SHexCer": "HexCer", "ASM": "SM", "LNAPE":"PE", "LNAPS":"PS", "NAPE":"PE", "NAPS":"PS",
                 "PE_Cer":"Cer",
                 "PI_Cer":"Cer", "ADGGA":"GL", "GM3":"GL", "DGDG":"GL", "DGGA":"GL", "MGDG":"GL", "SQDG":"GL",
                 "SHex":"SL", "STSE":"SL"}

def msdialnames_to_LC(df):
    fixednames = []
    for i, row in df.iterrows():
        name = row["Molecule Name"]
        if "|" in name:
            name = name.split("|")[1]
        name = fix_dash(name)
        name, _ = fix_name(name, "")
        fixednames.append(name)
    df["Molecule Name"] = fixednames
    return df


def fix_plasmalogens(name):
    if "TAG" in name:
        if "P_" in name:
            # Remove P- and add p to the end
            name = name.replace("P_", "")
            name = name + "p"
            # print(name)
        if "O_" in name:
            # Remove O- and add a to the end
            name = name.replace("O_", "")
            name = name + "a"
            # print(name)
        return name
    if "O_" in name:
        # Replace O_ with O-
        name = name.replace("O_", "O-")
        # After the first number after the first ":", add "a"
        if ":" in name:
            parts = name.split(":", 1)
            if len(parts) > 1:
                first_part = parts[0]
                second_part = parts[1]
                # Find the first number in the second part
                num = ""
                for char in second_part:
                    if char.isdigit():
                        num += char
                    else: break
                if num:
                    new_second_part = second_part.replace(num, num + "a", 1)
                    name = first_part + ":" + new_second_part
                    # Add the rest of the parts back in
        # print(name)
    if "P_" in name:
        # Replace P_ with O-
        name = name.replace("P_", "O-")
        # After the first number after the first ":", add "p"
        if ":" in name:
            parts = name.split(":", 1)
            if len(parts) > 1:
                first_part = parts[0]
                second_part = parts[1]
                # Find the first number in the second part
                num = ""
                for char in second_part:
                    if char.isdigit():
                        num += char
                    else: break
                if num:
                    new_second_part = second_part.replace(num, num + "p", 1)
                    name = first_part + ":" + new_second_part
                    # Add the rest of the parts back in
        # print(name)
    return name

def lipidlist_to_LC(df):
    outnames = []
    fixednames = []
    fixedclasses = []
    for i, row in df.iterrows():
        lclass = row["Molecule List Name"]
        name = row["Molecule Name"]
        name = fix_dash(name)

        name, lclass = fix_name(name, lclass)

        name = fix_plasmalogens(name)

        adduct = row["Precursor Adduct"]

        if adduct in adduct_translator[:, 0]:
            index = np.where(adduct_translator[:, 0] == adduct)[0][0]
            adduct_row = adduct_translator[index]
            adduct = adduct_row[-1]

        if adduct in adduct_translator[:, 1]:
            index = np.where(adduct_translator[:, 1] == adduct)[0][0]
            adduct_row = adduct_translator[index]
            adduct = adduct_row[-1]

        name, adduct = fix_d7(name, adduct)
        # print(name, adduct)
        fixednames.append(name)
        fixedclasses.append(lclass)

        outname = f"{name} {adduct}"
        outnames.append(outname)
    df["Molecule List Name"] = fixedclasses
    df["Molecule Name"] = fixednames
    df["Full Name Adduct"] = outnames
    return df


def fix_dash(name):
    if "-" in name:
        name = name.replace("-", "_")

    if "O3" in name:
        name = name.replace("O3", "3")
    if "O2" in name:
        name = name.replace("O2", "2")
    if "O4" in name:
        name = name.replace("O4", "4")
    if "O5" in name:
        name = name.replace("O5", "5")
    if "(2OH)" in name:
        name = name.replace("(2OH)", "")
    return name


def fix_name(name, lclass):
    if name.startswith("TG "):
        name = name.replace("TG", "TAG")
        lclass = "TAG"

    if name.startswith("DG "):
        name = name.replace("DG", "DAG")
        lclass = "DAG"

    if name.startswith("MG "):
        name = name.replace("MG", "MAG")
        lclass = "MAG"

    return name, lclass

def fix_adduct(adduct, incol, outcol):
    if adduct in adduct_translator[:, incol]:
        index = np.where(adduct_translator[:, incol] == adduct)[0][0]
        adduct_row = adduct_translator[index]
        adduct = adduct_row[outcol]
    return adduct


def fix_d7(name, adduct):
    if "d7" in name:
        adduct = adduct.replace("M", "M7H2")
        name= name.replace("(d7)", " + D7")

    if "d9" in name:
        adduct = adduct.replace("M", "M9H2")
        name = name.replace("(d9)", " + D9")

    return name, adduct


def count_tails(name, verbose=False):
    if "CoQ" in name:
        # Extract the number after CoQ
        tus = int(name.replace("CoQ", "").strip())
        tls = 5 * tus
        # print("CoQ detected", tls, tus)
        return "CoQ", tls, tus

    # Remove values of just a name or a name and then an adduct
    if " " not in name or ":" not in name:
        return "NA", -1, -1

    tls = 0
    tus = 0

    if ";FA" in name or "(FA" in name:
        # Extract the numbers after the FA and before )
        fa_matches = re.findall(r'FA (\d+:\d+)', name)
        for match in fa_matches:
            try:
                tl, tu = match.split(":")
                tl = int(tl)
                tu = int(tu)
                tls += tl
                tus += tu
            except Exception as e:
                if verbose:
                    print("Unable to parse FA match: ", match, name)
                continue
        if "(FA" in name:
            # Remove everything between ( and ) if it contains FA
            cleanname = re.sub(r'\(FA \d+:\d+\)', '', name)
        elif ";FA" in name:
            # Remove everything after ;FA if it contains FA
            cleanname = re.sub(r';FA \d+:\d+', '', name)
    else:
        cleanname = name

    if "(O-" in name:
        # Extract numbers after O- and before )
        o_matches = re.findall(r'O-(\d+:\d+)', name)
        for match in o_matches:
            try:
                tl, tu = match.split(":")
                tl = int(tl)
                tu = int(tu)
                tls += tl
                tus += tu
            except Exception as e:
                if verbose:
                    print("Unable to parse O- match: ", match, name)
                continue
        # Remove everything between ( and ) if it contains O-
        cleanname = re.sub(r'\(O-\d+:\d+\)', '', cleanname)

    try:
        tname = cleanname.split(" ")[1]
    except Exception as e:
        print("Unable to parse tails: ", name)
        raise e

    # fas = tname.split("_")
    fas = re.split('_|/', tname)

    # Remove '' from fas
    for f in fas:
        if f == '':
            continue
        if "O-" in f:
            f = f.replace("O-", "")
        if "P-" in f:
            f = f.replace("P-", "")
        if ";O2" in f:
            f = f.replace(";O2", "")
        if ";O3" in f:
            f = f.replace(";O3", "")
        if ";2O" in f:
            f = f.replace(";2O", "")
        if ";3O" in f:
            f = f.replace(";3O", "")
        if ";4O" in f:
            f = f.replace(";4O", "")
        if ";(2OH)" in f:
            f = f.replace(";(2OH)", "")
        if ";(3OH)" in f:
            f = f.replace(";(3OH)", "")
        if ";O" in f:
            f = f.replace(";O", "")
        if "N-" in f:
            f = f.replace("N-", "")
        if "-SN1" in f:
            f = f.replace("-SN1", "")
        if "-SN2" in f:
            f = f.replace("-SN2", "")
        if ";Hex" in f:
            f = f.replace(";Hex", "")
        if "(d7)" in f:
            f = f.replace("(d7)", "")
        if "(d9)" in f:
            f = f.replace("(d9)", "")

        tl = f.split(":")[0]  # Get length
        try:
            tu = f.split(":")[1]  # Get unsaturation
        except:
            if verbose:
                print("Unable to parse part 2: ", f, name, tname, fas)
            continue
            # exit()
        try:
            tl = int(tl)
        except:
            tls = -1
            if verbose:
                print("Unable to parse length 2: ", f, name, tname, fas)
            break
        try:
            tu = int(tu)
        except:
            tus = -1
            if verbose:
                print("Unable to parse unsaturation 2: ", f, name, tname, fas)
            break
        tls += tl
        tus += tu
    tl = tls
    tu = tus

    if tu == -1:
        if verbose:
            print("Unable to parse unsaturation: ", name, tname, fas)

    return tname, tl, tu


def set_tails(df, columnname="Molecule Name"):
    tails = []
    lengths = []
    unsats = []
    for i, row in df.iterrows():
        name = row[columnname]
        tname, tl, tu = count_tails(name, verbose=True)
        tails.append(tname)
        lengths.append(tl)
        unsats.append(tu)
    df["Tails"] = tails
    df["Tail Length"] = lengths
    df["Tail Unsaturation"] = unsats
    return df

def set_basic_tail_names(df, columnname="Molecule Name"):
    df = df.copy()
    tnames = []
    t1s = []
    t2s = []
    t3s = []
    t4s = []
    for i, row in df.iterrows():
        name = row[columnname]
        if "|" in name:
            name = name.split("|")[1]
        try:
            tname = name.split(" ", 1)[1]
        except Exception as e:
            tname = name
            print("Error parsing basic tail name:", name)
        t5 = ""
        if tname[0] == "(":
            # Extract stuff inside parentheses
            t5 = re.search(r'\((.*?)\)', tname).group(0)
            # Remove parentheses part from tname
            if t5:
                tname = tname.replace(t5, "")
                t5 = t5.replace("(", "").replace(")", "").strip()

        if "(FA" in tname:
            t5 = re.search(r'\(FA(.*?)\)', tname).group(0)
            if t5:
                tname = tname.replace(t5, "")
                t5 = t5.replace("(FA", "").replace(")", "").strip()

        fanames = re.split('_|/', tname)
        tnames.append(tname)
        t1 = fanames[0] if len(fanames) > 0 else ""
        t2 = fanames[1] if len(fanames) > 1 else ""
        t3 = fanames[2] if len(fanames) > 2 else ""
        t4 = fanames[3] if len(fanames) > 3 else ""

        if t5 != "":
            if t1 == "":
                t1 = t5
            elif t2 == "":
                t2 = t5
            elif t3 == "":
                t3 = t5
            elif t4 == "":
                t4 = t5

        t1s.append(t1)
        t2s.append(t2)
        t3s.append(t3)
        t4s.append(t4)

        # if t5 != "":
        #     print(t1, t2, t3, t4, t5, tname, name)
        #     exit()

    df["Tails"] = tnames
    df["T1"] = t1s
    df["T2"] = t2s
    df["T3"] = t3s
    df["T4"] = t4s

    return df

def get_fa_count(name):
    if "_" in name or "/" in name:
        return len(name.replace(" ", "").replace(")", "").replace("(", "").replace("/", "_").split("_"))
    else:
        return 1


def check_FA_count(df):
    df = df.copy()
    # set tails
    df = set_tails(df)
    # Iterate over df
    for i, row in df.iterrows():
        name = row["Molecule Name"]
        lclass = row["Molecule List Name"]
        if lclass in class_nfattyacids:
            nfa = class_nfattyacids[lclass]
            facount = get_fa_count(name)
            if facount > nfa:
                print("Too many FAs, dropping this row from df:", name, lclass)
                df = df.drop(i)
            if facount == 2 and nfa == 4:
                # Convert to 4FA format
                newname = convert_2FA_to_4FA(name, lclass)
                if newname is not None:
                    # Change name in df
                    df.at[i, "Molecule Name"] = newname

            if facount == 1 and nfa ==2:
                print("Warning: ", name, " in class ", lclass, " has 2 FAs, expected 1")
                # Drop it
                df = df.drop(i)
            # if facount != nfa:
            #     print(f"Warning: {name} in class {lclass} has {facount} FAs, expected {nfa}")
        # else:
        #     print(f"Warning: {lclass} not in class_nfattyacids")
    # Drop tails columns
    df = df.drop(columns=["Tails", "Tail Length", "Tail Unsaturation"])
    return df


def convert_2FA_to_4FA(name, lclass):
    tname, tl, tu = count_tails(name)
    if tl < 0 or tu < 0:
        print("Unable to convert, bad tail counts:", name)
        return None
    hg = name.split(" ")[0]
    combined_name = str(tl) + ":" + str(tu)
    newname = f"{hg} {combined_name}"
    print("Converting 2FA to SumFA for", name, newname)
    return newname


def remove_specific_problems(df):
    # Remove CL class with Na adduct
    # df = df[~((df["Molecule List Name"] == "CL") & (df["Full Name Adduct"].str.contains("Na")))]

    # Drop classes where Molecule List Name is one of these
    dropclasses = [#"CAR", "CE", "BA", "ST", 'FA', 'DAGDAG', 'SPB', 'BASulfate', 'BileAcid', 'HBMP', 'DHSph', 'SE','LNAPE',
         "Cer_AP", "Cer_HDS", 'Cer_NDS', 'Cer_NP', 'Cer_NS', 'FA'
        #'NAE', 'NAGly', 'NAGlySer', 'OxFA',
                   #'NAOrn', 'KLCAE','SSulfate', 'SL',
                   ]
    df = df[~df["Molecule List Name"].isin(dropclasses)]

    dropadducts = ["H2O", "Na"]#, "SMAGDAG", "DAGDAG", "BA", "VAE", ]
    # Drop if Full Name Adduct contains any of these as a substring
    for adduct in dropadducts:
        df = df[~df["Full Name Adduct"].str.contains(adduct)]

    dropfullnames = []  # "Cer 17:0;2/17:1;O [M+H]1+", 'Cer 20:0;2/26:0;O [M+H]1+', 'Cer 18:1;2/15:0 + D7 [M7H2+H]']
    df = df[~df["Full Name Adduct"].isin(dropfullnames)]

    classes = df['Molecule List Name'].unique()
    print("Remaining classes:", classes)

    # includeclasses = ["PC", "PE", "PS", "PI", "CL", "DAG", "TAG", "MAG", "SE", "PG", "LPG", "LPI", "LPC", "LCB", "LPE"]
    # includeclasses = ["SM"]
    #
    # # Which classes are not included
    # notincluded = set(classes) - set(includeclasses)
    # print("Not included classes:", notincluded)
    #
    # if len(includeclasses) > 0:
    #     df = df[df['Molecule List Name'].str.contains('|'.join(includeclasses), case=False, na=False)]

    return df


def map_ontology_from_lipid_name(lipid_name):
    """Map a lipid name to a lipid class label (e.g., PC, LPC, EtherPC)."""
    if lipid_name is None:
        return ""

    text = str(lipid_name).strip()
    if not text:
        return ""

    # Detect ether/plasmalogen notation before extracting the class token.
    is_ether = (
            text.startswith("O-")
            or text.startswith("P-")
            or text.startswith("O_")
            or text.startswith("P_")
            or text.startswith("ETHER")
            or " ETHER" in text
            or " P-" in text
            or " O-" in text
    )

    if "SM" in text and "FA" in text:
        return "ASM"

    # Keep primary annotation when names include pipe-delimited aliases.
    text = text.split(" ")[0].strip()

    # Output lipid class labels rather than broad ontology families.
    class_map = {
        'Cholesterol': 'ST',
        'CoQ1': 'CoQ',
        'CoQ2': 'CoQ',
        'CoQ3': 'CoQ',
        'CoQ4': 'CoQ',
        'CoQ5': 'CoQ',
        'CoQ6': 'CoQ',
        'CoQ7': 'CoQ',
        'CoQ8': 'CoQ',
        'CoQ9': 'CoQ',
        'CoQ10': 'CoQ',
        'CoQ11': 'CoQ',
        'CoQ12': 'CoQ',
        'CoQ13': 'CoQ',
        'PI-Cer': 'PI_Cer',
        'Cholesterol(d7)': 'ST',
        '25-hydroxycholecalciferol': 'Vitamin_D',
        'PE-Cer': 'PE_Cer',
        'Tocopherol': 'Vitamin_E',
    }

    if text in class_map:
        base_class = class_map[text]
    else:
        base_class = text

    if is_ether:
        return f"Ether{base_class}"

    return base_class

def apply_ontology_mapping(df, name_col="Metabolite name", output_col="Ontology"):
    df[output_col] = df[name_col].apply(map_ontology_from_lipid_name)
    return df

def make_sum_comp(name):
    if "_" or "/" in name:
        prefix = ""
        if "O-" in name:
            prefix = "O-"
        elif "P-" in name:
            prefix = "P-"
        specialpart = ""
        if ";" in name:
            specialpart = name.split(";")[1].strip()
            # Split off / or _ after the name
            specialpart = specialpart.split("/")[0].split("_")[0]
            specialpart = ";" + specialpart

        tname, tl, tu = count_tails(name, verbose=True)
        hname = name.split(" ")[0]
        newname = f"{hname} {prefix}{tl}:{tu}"
        name = newname + specialpart
    return name

def add_sump_comp(name):
    newname = make_sum_comp(name)
    return newname + "|" + name

def apply_sum_comp(df, namecol="Metabolite name", classcol="Ontology", adductcol="Adduct type"):
    # Apply to all
    # df["Sum Comp Name"] = df[namecol].apply(make_sum_comp)
    # # Create name as "Sum Comp Name" + " | " + original name
    # df = df.rename(columns={namecol: "Original Name"})
    # df[namecol] = df["Sum Comp Name"] + " | " + df["Original Name"]
    # df = df.drop(columns=["Sum Comp Name"])

    # Take only LPC, SM, and PC with [M+CH3COO]- adducts and apply make_sum_comp to them
    mask = df[adductcol] == "[M+CH3COO]-"
    mask2 = df[classcol].isin(["LPC", "SM", "PC"])
    if mask.any() and mask2.any():
        # Isolate LPC, SM, and PC [M+CH3COO]- adducts
        subdf = df[mask & mask2].copy()
        subdf[namecol] = subdf[namecol].apply(add_sump_comp)
        # Drop original LPC, SM, and PC [M+CH3COO]- rows
        df = df[~(mask & mask2)]
        # Append modified LPC, SM, and PC [M+CH3COO]- rows
        df = pd.concat([df, subdf], ignore_index=True)


    # mask = df[adductcol] == "[M+H]+"
    # mask2 = df[classcol] == "PC"
    # if mask.any() and mask2.any():
    #     # Isolate PC [M+H]+ adducts
    #     subdf = df[mask & mask2].copy()
    #     subdf[namecol] = subdf[namecol].apply(make_sum_comp)
    #     # Drop original PC [M+H]+ rows
    #     df = df[~(mask & mask2)]
    #     # Append modified PC [M+H]+ rows
    #     df = pd.concat([df, subdf], ignore_index=True)
    #
    # mask2 = df[classcol] == "SM"
    # if mask.any() and mask2.any():
    #     # Isolate SM [M+H]+ adducts
    #     subdf = df[mask & mask2].copy()
    #     subdf[namecol] = subdf[namecol].apply(make_sum_comp)
    #     # Drop original SM [M+H]+ rows
    #     df = df[~(mask & mask2)]
    #     # Append modified SM [M+H]+ rows
    #     df = pd.concat([df, subdf], ignore_index=True)
    return df

def get_charge(adduct):
    if adduct in ["[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M]+", "[M+H-H2O]+"]:
        return 1
    elif adduct in ["[M+2H]2+"]:
        return 2
    elif adduct in ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"]:
        return -1
    elif adduct in ["[M-2H]2-"]:
        return -2
    else:
        print("Unknown adduct", adduct)
        return 0

def add_product_ions(row, cutoff=0.1, spectrumkey="Ref Spec"):
    msmsspec = row[spectrumkey]
    # Remove final space if present
    if isinstance(msmsspec, str) and msmsspec.endswith(" "):
        msmsspec = msmsspec[:-1]
    # Parse spectrum
    msmsspec = [a.split(":") for a in msmsspec.split(" ")]
    try:
        msmsspec = np.array([(float(a[0]), float(a[1])) for a in msmsspec])
    except:
        print("Error parsing spectrum:", msmsspec)
        raise Exception("Error parsing spectrum")

    # Filter by cutoff
    b1 = msmsspec[:, 1] >= cutoff * np.max(msmsspec[:, 1])
    msmsspec = msmsspec[b1]

    mz = msmsspec[:, 0]
    newrows = []
    for m in mz:
        newrow = row.copy()
        charge = row["Precursor Charge"]
        if charge > 1:
            # Assume charge 1 for product ions
            charge = 1
        elif charge < -1:
            charge = -1
        newrow["Product Mz"] = m
        newrow["Product Charge"] = charge
        newrows.append(newrow)

    return newrows

def msdial_to_transitionlist(df, cutoff=0, spectrumkey="Ref Spec", rtwindow=4, remove_precursors=True):
    # Drop all where "Ontology" is Others or Unknown
    # df = df[~df["Ontology"].isin(["Others", "Unknown"])]
    newdf = df.copy()
    # Export Df with "Molecule List Name", "Molecule Name", "Precursor Adduct"
    newdf = newdf[["Ontology", "Metabolite name", "Adduct type", "Reference m/z", "Formula", "Average Rt(min)", spectrumkey]]
    newdf = newdf.rename(columns={"Ontology": "Molecule List Name", "Metabolite name": "Molecule Name",
                            "Adduct type": "Precursor Adduct", "Reference m/z": "Precursor Mz",
                                  "Average Rt(min)":"Explicit Retention Time", "Formula":"Molecule Formula"})

    charges = []
    for i, row in newdf.iterrows():
        prec_adduct = row["Precursor Adduct"]
        charge = get_charge(prec_adduct)
        charges.append(charge)
    newdf["Precursor Charge"] = charges

    # Set Explicit Retention Time Window to 4
    newdf["Explicit Retention Time Window"] = rtwindow

    # Sort by Molecule List Name and Molecule Name
    newdf = newdf.sort_values(by=["Molecule List Name", "Molecule Name"])

    newrows = []
    # Get product ion list
    for i, row in newdf.iterrows():
        prodrows = add_product_ions(row, cutoff=cutoff, spectrumkey=spectrumkey)
        newrows.extend(prodrows)

    # New rows to df
    proddf = pd.DataFrame(newrows)

    # Drop MS/MS spectrum column
    proddf = proddf.drop(columns=[spectrumkey], errors="ignore")
    proddf = pd.DataFrame(newrows)

    # If remove_precursors is True, drop rows where Product Mz is greater than 5 below Precursor Mz
    if remove_precursors:
        proddf = proddf[proddf["Product Mz"] < proddf["Precursor Mz"] - 5]

    return proddf

def cleanup_and_merge_msdial(posfile, negfile, posxml, negxml, rttol=0.2):
    # Load the data
    negdf = pd.read_csv(negfile)
    posdf = pd.read_csv(posfile)

    # Sort by Ontology and Name
    posdf = posdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)
    negdf = negdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)

    # Load the reference spectra from XML
    refnegdf = parse_xml(negxml)
    refposdf = parse_xml(posxml)

    # Add reference spectra to the dataframes
    negdf = add_ref_from_xml(negdf, refnegdf)
    posdf = add_ref_from_xml(posdf, refposdf)

    # Merge positive and negative dataframes
    posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

    # Create polarity columns and fill with "Positive" and "Negative"
    posdf["Polarity"] = "Positive"
    negdf["Polarity"] = "Negative"

    # Concatenate the dataframes
    mergeddf = pd.concat([posdf, negdf], ignore_index=True)

    return mergeddf, posdf, negdf

def cleanup_and_merge_mzmine(posfile, negfile, rttol=0.2):
    # Load the data
    negdf = pd.read_csv(negfile)
    posdf = pd.read_csv(posfile)

    # Sort by Ontology and Name
    posdf = posdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)
    negdf = negdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)

    # Merge positive and negative dataframes
    posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

    # Create polarity columns and fill with "Positive" and "Negative"
    posdf["Polarity"] = "Positive"
    negdf["Polarity"] = "Negative"

    # Concatenate the dataframes
    mergeddf = pd.concat([posdf, negdf], ignore_index=True)

    return mergeddf, posdf, negdf


def find_match_name(name, rt, refdf, rttol=0.1, verbose=False):
    newname = ""
    # Look for same name in negdf contained ahead of "|"
    rows = refdf[refdf["Metabolite name"].str.contains(name + "|", regex=False)]
    if len(rows) == 0:
        return newname
    if len(rows) == 1:
        refrt = rows.iloc[0]["Average Rt(min)"]
        refname = rows.iloc[0]["Metabolite name"]
        # If RTs are within rttol, assume same lipid
        if abs(rt - refrt) < rttol:
            newname = refname
    else:
        # Find one with the closest rt
        rtdiff = np.abs(rows["Average Rt(min)"] - rt)
        rows = rows[rtdiff == rtdiff.min()]
        if len(rows) > 1 and verbose:
            print("Warning Still multiple matches for", name, rt, rows[["Metabolite name", "Average Rt(min)"]])
        refrt = rows.iloc[0]["Average Rt(min)"]
        refname = rows.iloc[0]["Metabolite name"]
        if abs(rt - refrt) < rttol:
            newname = refname
    return newname

def merge_pos_neg(posdf, negdf, rttol=0.1):
    # Loop through posdf and look for names without "|"
    for i, row in posdf.iterrows():
        name = row["Metabolite name"]
        rtpos = row["Average Rt(min)"]
        if "|" not in name:
            newname = find_match_name(name, rtpos, negdf, rttol=rttol)
            if newname != "":
                posdf.at[i, "Metabolite name"] = newname
                print("Positive Updated name:", name, "to", newname)

    # Repeat in negative mode
    for i, row in negdf.iterrows():
        name = row["Metabolite name"]
        negrt = row["Average Rt(min)"]
        if "|" not in name:
            newname = find_match_name(name, negrt, posdf, rttol=rttol)
            if newname != "":
                negdf.at[i, "Metabolite name"] = newname
                print("Negative Updated name:", name, "to", newname)

    return posdf, negdf

def add_ref_spec_xml(df, xmlfile):
    refdf = parse_xml(xmlfile)
    print("Parsed", len(refdf), "spectra from XML file")
    df = add_ref_from_xml(df, refdf)
    return df

def add_ref_from_xml(df, refdf):
    # Match SpotID in refdf to Alignment ID in df
    for index, row in refdf.iterrows():
        spotid = row["SpotID"]
        spec = row["MS/MS spectrum ref"]
        # Find rows in df that match Alignment ID
        rows = df[df["Alignment ID"] == int(spotid)]
        if len(rows) == 0:
            continue
        for i, r in rows.iterrows():
            df.at[i, "Ref Spec"] = spec

    # Print how many rows have Ref Spec empty
    n_empty = len(df[df["Ref Spec"].isnull()])
    print("N rows without Ref Spec:", n_empty, "out of", len(df))
    return df

def parse_xml(xmlfile, exclude_precursor=True):
    tree = lxml.etree.parse(xmlfile)
    root = tree.getroot()
    # Extract spots
    spots = root.find("Spots")
    print("N Spots:", len(spots))

    rows = []
    # Loop through spots
    for spot in spots:
        row = {}
        spotid = -1
        adduct = ""
        name = ""
        prec_mz = -1
        spectrum = None
        specstring = ""
        for elem in spot:
            # Extract SpotID
            if elem.tag == "SpotId":
                spotid = elem.text
                # print("SpotID:", spotid)
                continue

            # If "Reference" element exists, strip
            if elem.tag == "Reference":
                adduct = elem.find("Adduct").text
                name = elem.find("Name").text
                prec_mz = elem.find("Mz").text
                # print("Reference:", adduct, name, mz)
                continue

            spectrum = elem.find("MatchedSpectrum")
            if spectrum is not None:
                # Find number of MatchSpectrumPeaks
                npeaks = len(spectrum.findall("MatchedSpectrumPeak"))
                if len(spectrum) == 0:
                    continue
                for subelem in spectrum:

                    # If Peak is "Reference", strip
                    if subelem.tag == "MatchedSpectrumPeak":
                        for subsubelem in subelem:
                            if subsubelem.tag == "Reference":
                                # Extract Mz and Intensity
                                mz = subsubelem.find("Mz").text
                                inten = subsubelem.find("Intensity").text

                                if exclude_precursor and npeaks > 1:
                                    if abs(float(mz) - float(prec_mz)) < 5:
                                        # print("Excluding precursor mz:", mz)
                                        continue
                                # Round mz to 4 decimal places
                                mz = round(float(mz), 4)
                                if str(mz) in specstring:
                                    continue
                                specstring += f"{mz}:{inten} "

        if spotid != -1 and specstring != "":
            row["SpotID"] = spotid
            row["Adduct"] = adduct
            row["Name"] = name
            row["Mz"] = prec_mz
            row["MS/MS spectrum ref"] = specstring
            rows.append(row)

    columns = ["SpotID", "Adduct", "Name", "Mz", "MS/MS spectrum ref"]
    df = pd.DataFrame(rows, columns=columns)
    return df

def clean_names(name):
    if "|" in name:
        name = name.split("|")[1]
    return name

def clean_df(df, column="Molecule name"):
    df = df.copy()
    cleanednames = []
    for i, row in df.iterrows():
        name = row[column]
        name = clean_names(name)
        cleanednames.append(name)
    df[column] = cleanednames
    return df

def encode_headgroups(df, classname="Ontology", rtcol_name="Average Rt(min)"):
    headgroups = df[classname].unique()
    if rtcol_name is None:
        sorted_hgs = sorted(headgroups)
    else:
        # Sort by average RT
        avg_rts = {}
        for hg in headgroups:
            subdf = df[df[classname] == hg]
            avg_rt = subdf[rtcol_name].mean()
            avg_rts[hg] = avg_rt
        # Sort headgroups by avg_rt
        sorted_hgs = sorted(avg_rts, key=avg_rts.get)
    # Create encoding dict
    encoding = {}
    for i, hg in enumerate(sorted_hgs):
        encoding[hg] = i + 1
    # print("Headgroup encoding:", encoding)
    # Add encoding column
    encodings = []
    for i, row in df.iterrows():
        hg = row[classname]
        enc = encoding[hg]
        encodings.append(enc)
    df["Headgroup Encoding"] = encodings

    # Create one-hot encoding columns
    for hg in headgroups:
        colname = f"HG_{hg}"
        df[colname] = df[classname].apply(lambda x: 1 if x == hg else 0)
    return df

def msdial_to_lipidlist(df):
    # Drop all where "Ontology" is Others or Unknown
    df = df[~df["Ontology"].isin(["Others", "Unknown"])]

    # Export Df with "Molecule List Name", "Molecule Name", "Precursor Adduct"
    df = df[["Ontology", "Metabolite name", "Adduct type"]]
    df = df.rename(columns={"Ontology": "Molecule List Name", "Metabolite name": "Molecule Name",
                            "Adduct type": "Precursor Adduct"})

    # Sort by Molecule List Name and Molecule Name
    df = df.sort_values(by=["Molecule List Name", "Molecule Name"])

    # Clean names
    df["Molecule Name"] = df["Molecule Name"].apply(clean_names)

    # Check FA count
    df = check_FA_count(df)

    # Convert to LC format
    df = lipidlist_to_LC(df)

    return df

def compress_similar_species(df, rttol=0.05):
    # Find all species that share the same name and adduct but have different retention times, compress them to the same retention time
    grouped = df.groupby(["Molecule Name"])
    rename_map = {}
    for (name, ), group in grouped:
        if len(group) > 1:
            # Check if all retention times are different and continue if they are all the same
            rts = group["Explicit Retention Time"].values
            if len(set(rts)) == 1:
                continue

            # If not, rename all with the same retention time with suffixes
            rt_groups = group.groupby("Explicit Retention Time")
            for rt, rt_group in rt_groups:
                if len(rt_group) >= 1:
                    # Select all from DF with the same name and rt within rttol
                    othermatchs = df[(df["Molecule Name"] == name) &
                                      (df["Explicit Retention Time"].between(rt - rttol, rt + rttol)) &
                                      (~df.index.isin(rt_group.index))]

                    otherrt = np.mean(othermatchs["Explicit Retention Time"].values) if len(othermatchs) > 0 else rt
                    avgrt = (rt + otherrt) / 2.0
                    # Round to 3 decimal places
                    avgrt = round(avgrt, 3)
                    for idx in rt_group.index:
                        if idx in rename_map:
                            continue
                        rename_map[idx] = avgrt

                    for om_idx in othermatchs.index:
                        if om_idx in rename_map:
                            continue
                        rename_map[om_idx] = avgrt

    df = df.copy()
    for idx, new_rt in rename_map.items():
        df.at[idx, "Explicit Retention Time"] = new_rt

    # Drop duplicate rows
    df = df.drop_duplicates(subset=["Molecule Name", "Precursor Adduct", "Explicit Retention Time", "Product Mz"])

    # Sort by Molecule Name, Precursor Adduct, Explicit Retention Time, Product Mz
    df = df.sort_values(by=["Molecule Name", "Precursor Adduct", "Explicit Retention Time", "Product Mz"])

    return df

def rename_overlapping_species(df, rttol=0.05):
    # Compress similar species first
    df = compress_similar_species(df, rttol=rttol)

    # Find all species that share the same name and adduct but have different retention times, rename them as _A _B, etc
    grouped = df.groupby(["Molecule Name"])
    rename_map = {}
    for (name,), group in grouped:
        if len(group) > 1:
            # Check if all retention times are different and continue if they are all the same
            rts = group["Explicit Retention Time"].values
            if len(set(rts)) == 1:
                continue

            # If not, rename all with the same retention time with suffixes
            rt_groups = group.groupby("Explicit Retention Time")
            suffixes = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            i=0
            for rt, rt_group in rt_groups:
                if len(rt_group) >= 1:
                    for j, idx in enumerate(rt_group.index):
                        if idx in rename_map:
                            continue
                        new_name = f"{name}_{suffixes[i]}"
                        rename_map[idx] = new_name

                    # Select all from DF with the same name and rt within rttol
                    othermatchs = df[(df["Molecule Name"] == name) &
                                      (df["Explicit Retention Time"].between(rt - rttol, rt + rttol)) &
                                      (~df.index.isin(rt_group.index))]
                    for om_idx in othermatchs.index:
                        if om_idx in rename_map:
                            continue
                        new_name = f"{name}_{suffixes[i]}"
                        rename_map[om_idx] = new_name

                    i+=1

    df = df.copy()
    for idx, new_name in rename_map.items():
        df.at[idx, "Molecule Name"] = new_name

    df = df.sort_values(by=["Molecule Name", "Precursor Adduct", "Explicit Retention Time", "Product Mz"])
    return df


def get_specstring(spectrum, remove_precursor=True, precursor_mz=None):
    spec = spectrum.peaks
    spec = np.array(spec)
    specstr = ""
    try:
        for mz, inten in zip(spec[:, 0], spec[:, 1]):
            mz = float(mz)
            if remove_precursor and precursor_mz is not None:
                if mz > precursor_mz - 5:
                    continue
            specstr += f"{mz:.4f}:{int(inten)} "
    except:
        print("Error in spectrum:", spec)
    specstr = specstr.strip()
    return specstr


def add_ref_spectra(df, libfile, rttol=0.1, remove_precursor=True):
    names = df["Metabolite name"].values
    redone_names = []
    for n in names:
        if "|" in n:
            n = n.split("|")[1]
        redone_names.append(n)
    df["Metabolite name Corr"] = redone_names
    from matchms.importing import load_from_msp
    lib = load_from_msp(libfile)
    outdata = []
    for l in lib:
        name = l.metadata["compound_name"]
        # print(name)
        if name in df["Metabolite name Corr"].values:
            adduct = l.metadata["adduct"]
            rt = l.metadata["retention_time"]
            libmz = l.metadata["precursor_mz"]

            # Find row in dfpos that matches name and adduct
            rows = df[df["Metabolite name Corr"] == name]
            rows = rows[rows["Adduct type"] == adduct]
            if len(rows) == 0:
                continue
            if len(rows) == 1:
                index = rows.index[0]
                rowrt = rows.at[index, "Average Rt(min)"]
                if np.abs(rowrt - rt) > rttol:
                    # print("RT mismatch for", name, adduct, rowrt, rt)
                    continue
                # else:
                #     print("Match for", name, adduct, rowrt, rt)

            # if len(rows) > 1:
            #     print("Multiple matches for", name, adduct, rows)

            outdata.append(l)
            specstr = get_specstring(l, remove_precursor=remove_precursor, precursor_mz=libmz)
            for index, row in rows.iterrows():
                # Check RT
                rowrt = row["Average Rt(min)"]
                if np.abs(rowrt - rt) > rttol:
                    continue
                #else:
                df.at[index, "Lib Spec"] = specstr

    print(len(outdata), len(df))
    # Count how many Lib Spec in df are empty
    nempty = np.sum(df["Lib Spec"].isna())
    print("N empty Lib Spec:", nempty)
    return df, outdata

def check_headgroup_fragments(df, fragdf, tol=0.01, intthresh=0.01):
    good = []
    goodnl = []
    for index, row in df.iterrows():
        class_name = row["Ontology"]
        adduct = row["Adduct type"]
        if "+" in adduct:
            ionmode = "positive"
        else:
            ionmode = "negative"

        fragments = fragdf[
            (fragdf["CompoundClass"] == class_name) &
            (fragdf["Adduct"] == adduct) &
            (fragdf["IonMode"] == ionmode)
        ]
        if fragments.empty:
            good.append(-1)
            goodnl.append(-1)
            # print("Not Found:", class_name, adduct, ionmode)
            continue

        msms_mz = row["MS/MS spectrum"]
        msms_mz = np.array([m.split(":") for m in msms_mz.split(" ")]).astype(float)
        b1 = msms_mz[:, 1] >= intthresh * np.max(msms_mz[:, 1])
        msms_mz = msms_mz[b1, :]

        prodmz = fragments["HeadgroupFragment_mz"].values[0]
        try:
            prodmz = np.array([prodmz]).astype(float)
        except:
            prodmz = np.array([float(f) for f in prodmz.split(",") if f])

        if prodmz is None or np.any(np.isnan(prodmz)) or len(prodmz) == 0:
            good.append(-1)
            # print("No Fragments:", class_name, adduct, ionmode)
        else:
            match_count = sum(1 for mz in msms_mz[:,0] if any(abs(mz - f) < tol for f in prodmz))
            good.append(match_count/len(prodmz))
            # if match_count == len(prodmz):
            #     good.append(1)
            # else:
            #     good.append(0)

        nlmz = fragments["HeadgroupNL_mz"].values[0]
        try:
            nlmz = np.array([nlmz]).astype(float)
        except:
            nlmz = np.array([float(f) for f in nlmz.split(",") if f])

        if nlmz is None or np.any(np.isnan(nlmz)) or len(nlmz) == 0:
            goodnl.append(-1)
            # print("No NL Fragments:", class_name, adduct, ionmode)
        else:
            nl_mz = row["Average Mz"] - msms_mz[:,0]
            match_count_nl = sum(1 for mz in nl_mz if any(abs(mz - f) < tol for f in nlmz))
            goodnl.append(match_count_nl/len(nlmz))
            # if match_count_nl == len(prodmz):
            #     goodnl.append(1)
            # else:
            #     goodnl.append(0)

        # print(f"Row {index}: Class={class_name} Adduct={adduct}, IonMode={ionmode}, Good={good[-1]} GoodNL={goodnl[-1]}")

    df["HG_Fragment_Match"] = np.array(good, dtype=float).round(2)
    df["HG_NL_Match"] = np.array(goodnl, dtype=float).round(2)

    # Set match counts > 1 to 1
    df.loc[df["HG_Fragment_Match"] > 1, "HG_Fragment_Match"] = 1.0
    df.loc[df["HG_NL_Match"] > 1, "HG_NL_Match"] = 1.0
    return df

def check_tail_fragments(df, taildf, tol=0.01, intthresh=0.01, allowedclasses=None):
    good = []
    goodnl = []
    df = set_basic_tail_names(df, columnname="Metabolite name")
    for index, row in df.iterrows():
        class_name = row["Ontology"]
        if allowedclasses is not None:
            if class_name not in allowedclasses:
                # Append -2 for all tails
                good.extend([-2, -2, -2, -2])
                goodnl.extend([-2, -2, -2, -2])
                continue
        adduct = row["Adduct type"]
        if "+" in adduct:
            ionmode = "positive"
        else:
            ionmode = "negative"

        t1 = row["T1"]
        t2 = row["T2"]
        t3 = row["T3"]
        t4 = row["T4"]
        tlist = [t1, t2, t3, t4]
        for tail in tlist:
            if pd.isna(tail) or tail == "":
                good.append(-2)
                goodnl.append(-2)
                continue
            fragments = taildf[
                (taildf["CompoundClass"] == class_name) &
                (taildf["Adduct"] == adduct) &
                (taildf["IonMode"] == ionmode) &
                (taildf["Tail"] == tail)]

            if fragments.empty:
                good.append(-1)
                goodnl.append(-1)
                # print("Not Found:", class_name, adduct, ionmode)
                continue

            msms_mz = row["MS/MS spectrum"]
            msms_mz = np.array([m.split(":") for m in msms_mz.split(" ")]).astype(float)
            b1 = msms_mz[:, 1] >= intthresh * np.max(msms_mz[:, 1])
            msms_mz = msms_mz[b1, :]

            # print(f"Row {index}, Tail {tail}: Class={class_name} Adduct={adduct}, IonMode={ionmode}")
            # print(fragments)

            prodmz = fragments["TailFragment_mz"].values[0]
            try:
                prodmz = np.array([prodmz]).astype(float)
            except:
                prodmz = np.array([float(f) for f in prodmz.split(",") if f])

            if prodmz is None or np.any(np.isnan(prodmz)) or len(prodmz) == 0:
                good.append(-1)
                # print("No Fragments:", class_name, adduct, ionmode)
            else:
                match_count = sum(1 for mz in msms_mz[:, 0] if any(abs(mz - f) < tol for f in prodmz))
                good.append(match_count / len(prodmz))


            nlmz = fragments["TailNL_mz"].values[0]
            try:
                nlmz = np.array([nlmz]).astype(float)
            except:
                nlmz = np.array([float(f) for f in nlmz.split(",") if f])

            if nlmz is None or np.any(np.isnan(nlmz)) or len(nlmz) == 0:
                goodnl.append(-1)
                # print("No NL Fragments:", class_name, adduct, ionmode)
            else:
                nl_mz = row["Average Mz"] - msms_mz[:, 0]
                match_count_nl = sum(1 for mz in nl_mz if any(abs(mz - f) < tol for f in nlmz))
                goodnl.append(match_count_nl / len(nlmz))

            print(
                f"Row {index}: Class={class_name} Adduct={adduct}, IonMode={ionmode}, Good={good[-1]} GoodNL={goodnl[-1]}")

    good = np.array(good, dtype=float).round(2)
    good[good > 1] = 1.0
    # Set all below 0 to NaN
    good[good < 0] = np.nan
    # Reshape good and goodnl to have one entry per row in df
    n_tails = 4
    n_rows = len(df)
    good = good.reshape((n_rows, n_tails))
    goodnl = np.array(goodnl, dtype=float).round(2)
    goodnl[goodnl > 1] = 1.0
    goodnl[goodnl < 0] = np.nan
    goodnl = goodnl.reshape((n_rows, n_tails))

    df["T1_Match"] = good[:, 0]
    df["T2_Match"] = good[:, 1]
    df["T3_Match"] = good[:, 2]
    df["T4_Match"] = good[:, 3]
    df["T1_NL_Match"] = goodnl[:, 0]
    df["T2_NL_Match"] = goodnl[:, 1]
    df["T3_NL_Match"] = goodnl[:, 2]
    df["T4_NL_Match"] = goodnl[:, 3]

    # Combine tail matches into overall Tail_Fragment_Match and Tail_NL_Match
    df["Tail_Fragment_Match"] = np.nansum(good, axis=1)
    df["Tail_NL_Match"] = np.nansum(goodnl, axis=1)

    return df

def categorize_nl(x):
    if pd.isna(x) or x == -1:
        return "No NL"
    if x == 1:
        return "All"
    if x == 0:
        return "None"
    return "Partial"

def merge_grades(grade1, grade2):
    # If grade2 is alphabetically first, swap
    if grade2 < grade1:
        grade1, grade2 = grade2, grade1
    a = grade1[0]
    b = grade2[0]
    c = "F"
    if a == "A" or b == "A":
        c = "A: Headgroup and All Tails Match"
    elif grade1 == "B: Headgroup and Some Tails Match" and grade2 == "B: Headgroup and Some Tails Match":
        c = "BB: Headgroup and Some Tails Match. Manually check for both."
    elif a == "B" and (grade2 == "C: Headgroup Match Only" or b=="D" or b=="F"):
        c = grade1
    elif a == "B" and grade2 == "C: All Tails Match Only":
        c = "A: Headgroup and All Tails Match"
    elif grade1 == "C: Headgroup Match Only" and grade2 == "C: Headgroup Match Only":
        c = grade1
    elif grade1 == "C: All Tails Match Only" and grade2 == "C: All Tails Match Only":
        c = grade1
    elif grade1 == "C: All Tails Match Only" and grade2 == "C: Headgroup Match Only":
        c = "A: Headgroup and All Tails Match"
    elif a == "C" and (b == "D" or b == "F"):
        c = grade1
    elif a == "D" and b == "F":
        c = grade1
    elif grade1 == grade2:
        c = grade1
    else:
        c = "Unable to Merge Grades: " + grade1 + " + " + grade2
    return c

def get_unique_names(df, namecol="Metabolite name", classcol="Ontology", rtcol="Average Rt(min)",
                     gradecol="Match_Grade", outlier_col="outlier", rttol=0.1):
    unique_names = {}
    for i, row in df.iterrows():
        name = row[namecol]
        classname = row[classcol]
        rt = row[rtcol]
        match_grade = row[gradecol]
        outlier = bool(row[outlier_col]) if outlier_col in df.columns else False

        if name not in unique_names:
            unique_names[name] = {"Class": classname, "RT": rt, "Match_Grade": match_grade, "Outlier": outlier}
        else:
            existing_rt = unique_names[name]["RT"]
            if abs(existing_rt - rt) > rttol:
                # Add a suffix to the name to make it unique
                new_name = f"{name}_rt{rt:.2f}"
                unique_names[new_name] = {"Class": classname, "RT": rt, "Match_Grade": match_grade, "Outlier": outlier}
            else:
                # Merge the grades
                existing_grade = unique_names[name]["Match_Grade"]
                newgrade = merge_grades(existing_grade, match_grade)
                unique_names[name]["Match_Grade"] = newgrade
                # If one is not an outlier, set to false
                if not outlier:
                    unique_names[name]["Outlier"] = False

    newdf = pd.DataFrame()
    for name in unique_names:
        classname = unique_names[name]["Class"]
        rt = unique_names[name]["RT"]
        newrow = {"Metabolite name": name, "Ontology": classname, "Average Rt(min)": rt,
                  "Match_Grade": unique_names[name]["Match_Grade"], "outlier": unique_names[name]["Outlier"]}
        newdf = pd.concat([newdf, pd.DataFrame([newrow])], ignore_index=True)
    return newdf


def detect_outliers_class(df, features, separate_adducts=False, contribs=True, sort=False):
    newdf = pd.DataFrame()
    for c in df['Ontology'].unique():
        # Get subset for class
        subdf = df[df['Ontology'] == c]

        # Clean up things that will cause issues with covariance matrix
        if len(subdf) < 3:
            newdf = pd.concat([newdf, subdf])
            continue

        if not separate_adducts:
            outputdf = detect_outliers_mahalanobis(subdf, features, contribs=contribs, sort=sort)
            newdf = pd.concat([newdf, outputdf])
        else:
            for adduct in subdf['Adduct type'].unique():
                subsubdf = subdf[subdf['Adduct type'] == adduct]
                if len(subsubdf) < 3:
                    newdf = pd.concat([newdf, subsubdf])
                    continue
                outputdf = detect_outliers_mahalanobis(subsubdf, features, contribs=contribs, sort=sort)
                newdf = pd.concat([newdf, outputdf])

    df = newdf.reset_index(drop=True)
    # Sort by outlier score
    df = df.sort_values(by='mahalanobis').reset_index(drop=True)

    # Print number of outliers detected
    print("Outlier detection complete. Number of outliers detected:", np.nansum(df['outlier']), "out of", len(df))
    return df


def mahalanobis_distance(data, robust=True, contrib="pca", dropcorr=None, checkrank=True):
    # Remove highly correlated variables
    to_drop = []
    if dropcorr is not None:
        df = pd.DataFrame(data)
        corr_matrix = df.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        to_drop = [column for column in upper.columns if any(upper[column] > dropcorr)]
        if len(to_drop) > 0:
            print("Dropping highly correlated variables:", to_drop)
        df = df.drop(columns=to_drop)
        data = df.values

    if checkrank:
        cov = np.cov(data, rowvar=False)
        is_full_rank = np.linalg.matrix_rank(cov) == cov.shape[0]
        # print(is_full_rank)
    else:
        is_full_rank = True

    if robust and is_full_rank:
        mcd = MinCovDet().fit(data)
        mean = mcd.location_
        cov = mcd.covariance_
    else:
        mean = np.mean(data, axis=0)
        cov = np.cov(data, rowvar=False)

    # try:
    inv_cov = np.linalg.pinv(cov)
    # except:
    #     inv_cov = np.linalg.pinv(cov)
    distances = np.array([mahalanobis(row, mean, inv_cov) for row in data])

    if contrib == "pca":
        # Calculate contributions using PCA method
        eigvals, V = np.linalg.eigh(cov)
        eigvals = np.maximum(eigvals, 1e-12)  # Prevent division by zero

        z = V.T @ (data - mean).T
        y = z.T / np.sqrt(eigvals)

        comp_contrib = np.square(y)
        contrib = comp_contrib @ (V.T **2)
        contrib = contrib / np.array(distances)[:, None]
    else:
        contrib = np.zeros_like(data)

    # Fill back in 0 for dropped columns
    if dropcorr is not None and len(to_drop) > 0:
        full_contrib = np.zeros((data.shape[0], data.shape[1] + len(to_drop)))
        keep_indices = [i for i in range(data.shape[1] + len(to_drop)) if i not in to_drop]
        for i, idx in enumerate(keep_indices):
            full_contrib[:, idx] = contrib[:, i]
        contrib = full_contrib

    return np.array(distances), contrib


def detect_outliers_mahalanobis(df, features, threshold=2, robust=True, verbose=False, contribs=True, sort=True, dropbadfeatures=True):
    if dropbadfeatures:
        # Drop any features with zero variance
        variances = df[features].var()
        valid_features = variances[variances > 0].index.tolist()
        # Sort features by variance
        valid_features = sorted(valid_features, key=lambda x: variances[x], reverse=True)

        # Check for two columns of valid features that are identical
        to_drop = set()
        for i in range(len(valid_features)):
            for j in range(i + 1, len(valid_features)):
                if np.all(df[valid_features[i]] == df[valid_features[j]]):
                    to_drop.add(valid_features[j])
        valid_features = [f for f in valid_features if f not in to_drop]

        # Remove class if not enough valid features
        if len(df) < len(valid_features) + 1:
            valid_features = valid_features[:len(df) - 1]
        features = valid_features

    # Calculate Mahalanobis distance for each point
    data = df[features].values
    outdf = df.copy()

    # Check for nan in data
    if np.any(np.isnan(data)):
        # Drop features with nan
        nan_cols = np.isnan(data).any(axis=0)
        nan_features = [features[i] for i in range(len(features)) if nan_cols[i]]
        print("Warning: NaN values found in data. Dropping features with NaN:", nan_features)
        data = df[[f for f in features if f not in nan_features]].values
        features = [f for f in features if f not in nan_features]

    md, contrib = mahalanobis_distance(data, robust=robust)
    outdf['mahalanobis'] = md
    # Print if any nan in md
    if np.any(np.isnan(md)):
        print("Warning: NaN values found in Mahalanobis distances. These rows will be marked as outliers.")
        # Set NaN distances to large value
        outdf.loc[np.isnan(md), 'mahalanobis'] = -10

    # Calculate mean and std dev of distances
    pmd = md[md<np.percentile(md, 95)]
    if len(pmd) > 0:
        md = pmd
    mean_md = np.mean(md)  # Exclude extreme outliers for mean/std calculation
    std_md = np.std(md)
    # Define outlier threshold
    threshold = mean_md + threshold * std_md

    outdf['outlier'] = outdf['mahalanobis'] > threshold

    if contribs:
        # Add contribution columns to outdf with feature_contrib name
        for i, feature in enumerate(features):
            outdf[f'{feature}_contrib'] = contrib[:, i]
        # Max contribution
        outdf["Top contrib"] = np.argmax(contrib, axis=1)
        outdf["Top contrib feature"] = outdf["Top contrib"].apply(lambda x: features[x])
        # Drop top contrib column
        outdf = outdf.drop(columns=["Top contrib"])

    if sort:
        # Sort by outlier score
        outdf = outdf.sort_values(by='mahalanobis').reset_index(drop=True)

    if verbose and sum(outdf['outlier']) != 0:
        print("Mahalanobis outlier detection complete. Number of outliers detected:", sum(outdf['outlier']), "of", len(outdf))
        testdf = outdf[outdf['outlier']][["Metabolite name"] + features + [f"{f}_contrib" for f in features]]
        print(testdf.to_string())
    return outdf


def class_network_analysis(df, min_count=3, header="Class_", add_all_cols=False, tol=0.05):
    graphs = []
    outdf = pd.DataFrame()
    for c in df['Ontology'].unique():
        # Get subset for class
        subdf = df[df['Ontology'] == c].reset_index(drop=True)

        if len(subdf) < min_count:
            outdf = pd.concat([outdf, subdf])
            continue

        newdf, G = network_analysis(subdf, header=header, add_all_cols=add_all_cols, mz_tol=tol)
        G.graph["Title"]=c
        graphs.append(G)
        outdf = pd.concat([outdf, newdf])
    outdf = outdf.reset_index(drop=True)
    return outdf, graphs

def correct_tail_network_columns(row, header):
    # Find whether the header is T1, T2, T3, or T4
    tails = ["T1", "T2", "T3", "T4"]
    tailmatches = []
    for i, tail in enumerate(tails):
        if row[tail] == header:
            tailmatches.append(tail)

    for t in tailmatches:
        # Find all columns that start with header
        tail_cols = [col for col in row.index if col.startswith(header)]

        # Find columns that start with header and replace with tail
        for col in tail_cols:
            if col.startswith(header):
                newcol = col.replace(header, t)
                row[newcol] = row[col]
                row = row.drop(labels=[col])

    return row

def merge_tail_rows(df):
    # Find all rows with same Alignment ID
    outdf = pd.DataFrame()
    grouped = df.groupby("Alignment ID")
    for name, group in grouped:
        if len(group) == 1:
            outdf = pd.concat([outdf, group])
            continue
        # Merge rows
        newrow = group.iloc[0].copy()
        for col in group.columns:
            if col in ["Alignment ID", "Metabolite name"]:
                continue
            # If any value in group[col] is not NaN, take that value
            values = group[col].dropna().unique()
            if len(values) > 0:
                newrow[col] = values[0]
        outdf = pd.concat([outdf, pd.DataFrame([newrow])])
    outdf = outdf.reset_index(drop=True)
    return outdf


def tail_network_analysis(df, mz_tol=0.01, columnname="Metabolite name", min_count=3, header="Tail_", add_all_cols=False):
    graphs = []
    outdf = pd.DataFrame()
    if "T1" not in df.columns:
        print("Creating Tail Names")
        df = set_basic_tail_names(df, columnname=columnname)
    tail_cols = ["T1", "T2", "T3", "T4"]
    uniquetails = set()
    for col in tail_cols:
        tails = df[col].unique()
        for t in tails:
            if pd.isna(t) or t == "":
                continue
            uniquetails.add(t)
    # Sort by tail name
    uniquetails = sorted(uniquetails)
    for tail in uniquetails:
        # Get subset for tail
        subdf = df[(df["T1"] == tail) | (df["T2"] == tail) | (df["T3"] == tail) | (df["T4"] == tail)].reset_index(drop=True)

        if len(subdf) < min_count:
            outdf = pd.concat([outdf, subdf])
            continue

        newdf, G = network_analysis(subdf, mz_tol=mz_tol, header=tail+"_", add_all_cols=add_all_cols)
        G.graph["Title"]=tail
        graphs.append(G)

        # Apply correct_tail_network_columns to each row in newdf
        newdf = newdf.apply(lambda row: correct_tail_network_columns(row, tail), axis=1)

        outdf = pd.concat([outdf, newdf])

    # For every unique Alignment ID
    outdf = merge_tail_rows(outdf)

    # Create an anomly score by summing tail network scores
    tail_scores = []
    tail_zscores = []
    for index, row in outdf.iterrows():
        ascore = 0.0
        zscore = 0.0
        tails = ["T1", "T2", "T3", "T4"]
        for tail in tails:
            score_col = f"{tail}_Anomaly_Score"
            zscore_col = f"{tail}_z_Anomaly_Score"
            if score_col in outdf.columns:
                val = row[score_col]
                if not pd.isna(val):
                    ascore += val
            if zscore_col in outdf.columns:
                zval = row[zscore_col]
                if not pd.isna(zval):
                    zscore += zval
        tail_scores.append(ascore)
        tail_zscores.append(zscore)

    outdf[header + "Anomaly_Score"] = tail_zscores
    outdf[header + "z_Anomaly_Score"] = robust_z_list(tail_zscores)

    return outdf, graphs



def robust_z(values_dict):
    # returns dict of z-scores
    vals = np.array(list(values_dict.values()), dtype=float)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    eps = 1e-9
    z = 0.6745 * (vals - med) / (mad + eps)
    return dict(zip(values_dict.keys(), z))

def robust_z_list(values):
    vals = np.array(values, dtype=float)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    eps = 1e-9
    z = 0.6745 * (vals - med) / (mad + eps)
    return z

def create_msms_score_matrix(df, mz_tol=0.01):
    # Create a score matrix based on MS/MS similarity
    spectra = []
    nlspectra = []
    for i, row in df.iterrows():
        if pd.isna(row['MS/MS spectrum']):
            raise ValueError(f"Row {i} has no MS/MS spectrum data.")
        text = row['MS/MS spectrum']
        data = text.split(" ")
        data2 = np.array([d.split(":") for d in data]).astype(float)
        mz = row['Average Mz']
        metadata = {"precursor_mz": mz}
        nldata = data2.copy()
        nldata[:, 0] = mz - data2[:, 0]
        # Sort nldata by m/z
        nldata = nldata[np.argsort(nldata[:, 0])]
        # Take only nldata > 10
        nldata = nldata[nldata[:, 0] > 10]
        spec = mms.Spectrum(data2[:, 0], data2[:, 1], metadata=metadata)
        spectra.append(spec)
        nlspec = mms.Spectrum(nldata[:, 0], nldata[:, 1], metadata=metadata)
        nlspectra.append(nlspec)

    score_matrix = mms.calculate_scores(references=spectra, queries=spectra,
                                        similarity_function=mms.similarity.CosineGreedy(),
                                        is_symmetric=True, ).scores.to_array()
    score_matrix = np.array(score_matrix.tolist())[:, :, 0]

    nlscore_matrix = mms.calculate_scores(references=nlspectra, queries=nlspectra,
                                          similarity_function=mms.similarity.CosineGreedy(),
                                          is_symmetric=True).scores.to_array()
    nlscore_matrix = np.array(nlscore_matrix.tolist())[:, :, 0]

    score_matrix += nlscore_matrix
    score_matrix /= 2.0
    return score_matrix

def network_analysis(df, mz_tol=0.01, title="", header="", add_all_cols=False):
    # Create score matrix
    score_matrix = create_msms_score_matrix(df, mz_tol=mz_tol)

    # Create graph from score matrix
    G = nx.Graph()
    G.add_weighted_edges_from([ (i, j, score_matrix[i,j]) for i in range(score_matrix.shape[0]) for j in range(i, score_matrix.shape[1]) if score_matrix[i,j] > 0.1 ])

    # Add node labels
    labels = {i: str(i) + "~"+ df.iloc[i]['Metabolite name'] for i in range(len(df))}
    nx.relabel_nodes(G, labels, copy=False)

    # Drop self-loops
    G.remove_edges_from(nx.selfloop_edges(G))

    # Add distance attribute for shortest-path calculations
    for u, v, data in G.edges(data=True):
        data['distance'] = 1.0 - data['weight']

    # Measure a variety of features of the graph
    strength = dict(G.degree(weight='weight'))  # s_i
    degree = dict(G.degree())  # unweighted degree
    # clust_w = nx.clustering(G, weight='weight')  # weighted clustering coeff
    # betw = nx.betweenness_centrality(G, weight='distance')  # uses shortest-path distances
    # close = nx.closeness_centrality(G, distance='distance')  # distance-based closeness
    # pr = nx.pagerank(G, weight='weight')  # weighted PageRank
    tri = nx.triangles(G)  # triangle count per node
    # ev = nx.eigenvector_centrality(G, weight='weight', tol=1e-04)  # eigenvector centrality
    # ----- SAFETY GUARD FOR EMPTY / EDGELESS GRAPHS -----
    # if G.number_of_nodes() < 2 or G.number_of_edges() == 0:
    #     ev = dict.fromkeys(G.nodes(), 0.0)
    # else:
    #     try:
    #         ev = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
    #     except nx.PowerIterationFailedConvergence:
    #         ev = dict.fromkeys(G.nodes(), 0.0)

    # Network features dict
    n_features = {
        'strength': strength,
        'degree': degree,
        # 'clust_w': clust_w,
        # 'betw': betw,
        # 'close': close,
        # 'pr': pr,
        'tri': tri,
        # 'ev': ev,
    }

    z_features = {
        'z_strength': robust_z(strength),
        'z_degree': robust_z(degree),
        # 'z_clust_w': robust_z(clust_w),
        # 'z_betw': robust_z(betw),
        # 'z_close': robust_z(close),
        # 'z_pr': robust_z(pr),
        'z_tri': robust_z(tri),
        # 'z_ev': robust_z(ev),
    }

    # Simple aggregate anomaly score A_i = sum |z_ik|
    nodes = list(G.nodes())
    A = {}
    for n in nodes:
        A[n] = sum(zf[n] for zf in z_features.values())
    # Add anomaly score as z_features
    z_features['Anomaly_Score'] = A
    z_features['z_Anomaly_Score'] = robust_z(A)
    outdf = df.copy()

    # Set these as columns in df
    for i, row in df.iterrows():
        name = str(i) + "~" + row['Metabolite name']
        if add_all_cols:
            # Apply to outdf
            for feat_name, feat_dict in n_features.items():
                outdf.at[i, header+feat_name] = float(feat_dict.get(name, 0.0))
            for zfeat_name, zfeat_dict in z_features.items():
                outdf.at[i, header+zfeat_name] = float(zfeat_dict.get(name, 0.0))
        else:
            # Only add anomaly scores
            outdf.at[i, header+'Anomaly_Score'] = float(z_features['Anomaly_Score'].get(name, 0.0))
            outdf.at[i, header+'z_Anomaly_Score'] = float(z_features['z_Anomaly_Score'].get(name, 0.0))

        if name in G:
            # Apply to nodes in G
            attrs = {}
            for feat_name, feat_dict in n_features.items():
                attrs[feat_name] = float(feat_dict.get(name, 0.0))
            for zfeat_name, zfeat_dict in z_features.items():
                attrs[zfeat_name] = float(zfeat_dict.get(name, 0.0))
            attrs['Adduct type'] = row['Adduct type']
            attrs['Class'] = row['Ontology']
            G.nodes[name].update(attrs)
    G.graph['Title'] = title
    return outdf, G

def get_fingerprint(mols, fingerprinttype="morgan"):
    mols = [Chem.RemoveAllHs(mol) for mol in mols]
    if fingerprinttype == "morgan":
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024, useChirality=False, useFeatures=True) for mol in mols]
    elif fingerprinttype == "rdkit":
        fps = [Chem.RDKFingerprint(mol) for mol in mols]
    elif fingerprinttype == "sparse":
        fpgen = AllChem.GetRDKitFPGenerator()
        fps = [fpgen.GetSparseFingerprint(mol) for mol in mols]
    elif fingerprinttype == "pattern":
        fps = [Chem.PatternFingerprint(mol) for mol in mols]
    elif fingerprinttype == "atompair":
        fpgen = AllChem.GetAtomPairGenerator()
        fps = [fpgen.GetSparseFingerprint(mol) for mol in mols]
    elif fingerprinttype == "torsion":
        fpgen = AllChem.GetTopologicalTorsionGenerator()
        fps = [fpgen.GetSparseFingerprint(mol) for mol in mols]
    elif fingerprinttype == "maccs":
        fps = [MACCSkeys.GenMACCSKeys(mol) for mol in mols]
    else:
        raise ValueError(f"Unknown fingerprint type: {fingerprinttype}")
    return fps

def get_similarity(fp1, fp2, simtype="tanimoto"):
    if simtype == "tanimoto":
        sim = DataStructs.TanimotoSimilarity(fp1, fp2)
    elif simtype == "dice":
        sim = DataStructs.DiceSimilarity(fp1, fp2)
    elif simtype == "cosine":
        sim = DataStructs.CosineSimilarity(fp1, fp2)
    elif simtype == "kulczynski":
        sim = DataStructs.KulczynskiSimilarity(fp1, fp2)
    elif simtype == "mcconnaughey":
        sim = DataStructs.McConnaugheySimilarity(fp1, fp2)
    elif simtype == "rogotgoldberg":
        sim = DataStructs.RogotGoldbergSimilarity(fp1, fp2)
    else:
        raise ValueError(f"Unknown similarity type: {simtype}")
    return sim

def gen_smiles_similarity_matrix(smiles_list, smiles_list2=None, simtype="tanimoto", fingerprinttype="morgan"):
    # Convert to mols
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]

    # Convert to fingerprints
    fps = get_fingerprint(mols, fingerprinttype=fingerprinttype)

    # Get similarity matrix
    if smiles_list2 is None:
        similarity_matrix = np.zeros((len(fps), len(fps)))
        for i in range(len(fps)):
            similarity_matrix[i, i] = 1.0
            for j in range(i + 1, len(fps)):
                sim = get_similarity(fps[i], fps[j], simtype=simtype)
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
    else:
        mols2 = [Chem.MolFromSmiles(s) for s in smiles_list2]
        # Loop over mol and set isotope to 0
        for mol in mols:
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            # AllChem.MMFFOptimizeMolecule(mol)
        for mol in mols2:
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            # AllChem.MMFFOptimizeMolecule(mol)


        fps2 = get_fingerprint(mols2, fingerprinttype=fingerprinttype)
        similarity_matrix = np.zeros((len(fps), len(fps2)))
        for i in range(len(fps)):
            for j in range(len(fps2)):
                sim = get_similarity(fps[i], fps2[j], simtype=simtype)
                similarity_matrix[i, j] = sim
    return similarity_matrix


def smiles_analysis(df, plot=False, cutoff=0.7, simtype="tanimoto", fingerprinttype="morgan"):
    outdf = df.copy()

    starttime = time.perf_counter()
    # Create network graph from SMILES
    smiles = df["SMILES"].dropna().unique()
    similarity_matrix = gen_smiles_similarity_matrix(smiles, simtype=simtype, fingerprinttype=fingerprinttype)

    print(f"Calculated similarity matrix, Time taken: {time.perf_counter() - starttime:.2f} seconds")
    # Convert to network graph using a similarity threshold of 0.1
    G = nx.Graph()
    G.add_weighted_edges_from(
        [(i, j, similarity_matrix[i, j]) for i in range(similarity_matrix.shape[0]) for j in
         range(i, similarity_matrix.shape[1]) if
         similarity_matrix[i, j] > cutoff])

    # Add node labels
    # Loop through df and find the index of each SMILES in the original df, then use that index to get the Metabolite name and add it to the label
    labels = {i: smiles[i] for i in range(len(smiles))}
    nx.relabel_nodes(G, labels, copy=False)

    # Drop self-loops
    G.remove_edges_from(nx.selfloop_edges(G))

    # Measure a strength for each in the graph
    strength = dict(G.degree(weight='weight'))

    mapping = {}
    # Set these as columns in df
    for i, row in df.iterrows():
        sval = row['SMILES']
        if pd.isna(sval):
            continue
        mapping[sval] = "0~"+row['Metabolite name']

        # Apply to outdf
        outdf.at[i, "SMILES_strength"] = float(strength.get(sval, 0.0))

        if sval in G:
            # Apply to nodes in G
            attrs = {}
            attrs['Strength'] = strength.get(sval, 0.0)
            attrs['Adduct type'] = row['Adduct type']
            attrs['Class'] = row['Ontology']
            attrs['Name'] = row['Metabolite name']
            G.nodes[sval].update(attrs)
    # G.graph['Title'] = title

    # Adjust node names to be the Metabolite name instead of the SMILES
    nx.relabel_nodes(G, mapping, copy=False)

    print(f"Created SMILES network graph, Time taken: {time.perf_counter() - starttime:.2f} seconds")

    if plot:
        # Plot the network graph
        network_plot(G, colortarget="Class")
        plt.show()
    return outdf

def network_plot(G, axis=None, colortarget="z_Anomaly_Score"):
    if axis is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        ax = axis

    ax.set_title(G.graph.get("Title", "Network Plot"))
    # Set up position layout
    pos = nx.spring_layout(G, k=0.5, iterations=50)

    # edges on the specific Axes
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    widths = [max(0.1, float(np.round(w, 2))) for w in weights]
    nx.draw_networkx_edges(G, pos, width=widths, ax=ax)

    # nodes as PathCollection on the specific Axes
    node_list = list(G.nodes())
    nc = nx.draw_networkx_nodes(G, pos, nodelist=node_list, ax=ax)

    # Set color as colortarget attribute
    colorvals = np.array([G.nodes[n].get(colortarget, 0.0) for n in node_list])
    # If colorvals are not floats, encode them as floats
    if not np.issubdtype(colorvals.dtype, np.floating):
        unique_vals = list(set(colorvals))
        # Sort alphanumerically
        unique_vals = sorted(unique_vals, key=lambda x: (str(type(x)), str(x)))
        val_to_float = {v: i for i, v in enumerate(unique_vals)}
        colorvals = np.array([val_to_float[v] for v in colorvals], dtype=float)
        # prepare label mapping
        labels_text = {node: (
            f"{node.split('~')[1]}\n"
            f"{colortarget}: {G.nodes[node].get(colortarget, 'N/A')}\n"
            f"Adduct: {G.nodes[node].get('Adduct type', 'N/A')}\n"
        ) for node in node_list}
        cmap = mpl.cm.tab20
    else:
        cmap = mpl.cm.viridis
        # prepare label mapping
        try:
            labels_text = {node: (
                f"{node.split('~')[1]}\n"
                f"{colortarget}: {G.nodes[node].get(colortarget, 0.0):.2f}\n"
                f"Adduct: {G.nodes[node].get('Adduct type', 'N/A')}\n"
            ) for node in node_list}
        except:
            labels_text = {node: str(i) for i, node in enumerate(node_list)}

    norm = mpl.colors.Normalize(vmin=np.amin(colorvals), vmax=np.amax(colorvals))

    node_colors = cmap(norm(colorvals))
    nc.set_facecolor(node_colors)

    # labels on the same Axes
    nx.draw_networkx_labels(G, pos, labels={n: n.split("~")[1] for n in node_list}, font_size=8, ax=ax)

    # attach cursor to this collection
    cursor = mplcursors.cursor(nc, hover=True)

    @cursor.connect("add")
    def on_add(sel, node_list=node_list, labels_text=labels_text, nc=nc):
        # robustly extract a single integer index
        idx = getattr(sel, "index", None)
        if idx is None:
            idx = getattr(getattr(sel, "target", None), "index", None)

        # if index is array-like, take first valid element
        if isinstance(idx, (np.ndarray, list, tuple)):
            if len(idx) == 0:
                idx = None
            else:
                try:
                    idx = int(np.atleast_1d(idx)[0])
                except Exception:
                    idx = None

        # fallback: nearest point by mouse coordinates
        if idx is None:
            try:
                offsets = nc.get_offsets()
                mx, my = sel.mouseevent.xdata, sel.mouseevent.ydata
                if mx is None or my is None:
                    return
                d = np.sum((offsets - np.array([mx, my])) ** 2, axis=1)
                idx = int(np.argmin(d))
            except Exception:
                return

        # validate index bounds
        if idx < 0 or idx >= len(node_list):
            return

        node = node_list[int(idx)]
        sel.annotation.set_text(labels_text[node])


def make_lipid_summary_plots(posdf, negdf):

    plt.figure(figsize=(8, 6))
    # Make plots of Average RT versus reference RT
    plt.subplot(221)
    x1 = posdf["Average Rt(min)"]
    y1 = posdf["Reference RT"]

    x2 = negdf["Average Rt(min)"]
    y2 = negdf["Reference RT"]

    plt.scatter(x1, y1, label="Positive", alpha=0.5, color="blue")
    plt.scatter(x2, y2, label="Negative", alpha=0.5, color="red")
    plt.plot([0, max(max(x1), max(x2))], [0, max(max(x1), max(x2))], color="black", linestyle="--")
    plt.xlabel("Average RT (min)")
    plt.ylabel("Reference RT (min)")
    plt.title("Average RT vs Reference RT")
    plt.legend()
    plt.grid()

    # Make violin plots of RT diff
    plt.subplot(222)
    rt_diff_pos = posdf["Average Rt(min)"] - posdf["Reference RT"]
    rt_diff_neg = negdf["Average Rt(min)"] - negdf["Reference RT"]

    plt.violinplot([rt_diff_pos.dropna(), rt_diff_neg.dropna()], positions=[1, 2], showmeans=True)
    plt.xticks([1, 2], ["Positive", "Negative"])
    plt.ylabel("RT Difference (min)")
    plt.title("RT Difference Distribution")
    plt.grid()

    # Make violin plots of Mz diff versus reference Mz
    plt.subplot(223)
    mz_diff_pos = posdf["Average Mz"] - posdf["Reference m/z"]
    mz_diff_neg = negdf["Average Mz"] - negdf["Reference m/z"]

    plt.violinplot([mz_diff_pos.dropna(), mz_diff_neg.dropna()], positions=[1, 2], showmeans=True)
    plt.xticks([1, 2], ["Positive", "Negative"])
    plt.ylabel("Mz Difference")
    plt.title("Mz Difference Distribution")
    plt.grid()

    # Violin plots of Simple dot product, Weighted dot product, Reverse dot product, and Total Score
    plt.subplot(224)
    sdp_pos = posdf["Simple dot product"]
    sdp_neg = negdf["Simple dot product"]
    wdp_pos = posdf["Weighted dot product"]
    wdp_neg = negdf["Weighted dot product"]
    rdp_pos = posdf["Reverse dot product"]
    rdp_neg = negdf["Reverse dot product"]
    ts_pos = posdf["Total score"]
    ts_neg = negdf["Total score"]

    plt.violinplot([sdp_pos.dropna(), sdp_neg.dropna(), wdp_pos.dropna(), wdp_neg.dropna(),
                    rdp_pos.dropna(), rdp_neg.dropna(), ts_pos.dropna(), ts_neg.dropna()],
                   positions=[1, 2, 3, 4, 5, 6, 7, 8], showmeans=True)
    plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["SDP Pos", "SDP Neg", "WDP Pos", "WDP Neg",
                                        "RDP Pos", "RDP Neg", "TS Pos", "TS Neg"], rotation=45)
    plt.ylabel("Score")
    plt.title("MS/MS Score Distributions")
    plt.grid()

    plt.tight_layout()
    plt.show()

def fix_inhomogeneous_rt(rtdata):
    lengths = [len(rt[0]) for rt in rtdata]
    # If all lengths are the same, return as is
    if len(set(lengths)) == 1:
        return np.array(rtdata)
    # Otherwise, interpolate to min to max range
    min_time = min([min(rt[0]) for rt in rtdata])
    max_time = max([max(rt[0]) for rt in rtdata])
    spacing = rtdata[0][0][1] - rtdata[0][0][0]
    new_times = np.arange(min_time, max_time, spacing)
    new_rtdata = []
    for times, intensities in rtdata:
        new_intensities = np.interp(new_times, times, intensities, left=0, right=0)
        new_rtdata.append((new_times, new_intensities))
    return np.array(new_rtdata)

def get_rt(rtrows, tcol="Times", icol="Intensities", mode="apex"):
    rtdata = []
    for index, row in rtrows.iterrows():
        times = [float(x) for x in row[tcol].split(",")]
        intensities = [float(x) for x in row[icol].split(",")]
        rtdata.append((times, intensities))
    rtdata = fix_inhomogeneous_rt(rtdata)
    times = rtdata[0,0]
    summed_eic = np.sum(rtdata[:,1], axis=0)
    rtsummed = np.transpose([times, summed_eic])

    com_rt = center_of_mass(rtsummed, relative_cutoff=0.1)
    if mode == "apex":
        max_index = np.argmax(rtsummed[:,1])
        max_pos = rtsummed[max_index,0]
        com_rt = (max_pos, com_rt[1])

    return float(com_rt[0]), float(com_rt[1])

def update_rt(df, cdf, min_window=3.0, max_window=5.0, mode="apex", nstds=8, forcedrts=None):
    # Find all unique molecule names in transition list
    unique_names = df["Molecule Name"].unique()
    outdf = df.copy()
    # Set Explicity Retention Time and Window as floats to avoid dtype issues
    outdf["Explicit Retention Time"] = outdf["Explicit Retention Time"].astype(float)
    outdf["Explicit Retention Time Window"] = outdf["Explicit Retention Time Window"].astype(float)

    for name in unique_names:
        if forcedrts is not None and name in forcedrts:
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time"] = forcedrts[name]
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time Window"] = min_window
            print(f"Manually setting {name} to RT: {forcedrts[name]:.2f} with window: {min_window:.2f}")
            continue

        # name = row["Molecule Name"]
        rt_row = cdf[cdf["PeptideModifiedSequence"] == name]
        # adduct = row["Precursor Adduct"][:-1]  # remove + or - at end
        # rt_row = rt_row[rt_row["PrecursorCharge"] == adduct]
        if not rt_row.empty:
            # df.at[index, "Explicit Retention Time"] = rt_row["Retention Time (min)"].values[0]
            # df.at[index, "Explicit Retention Time Window"] = rt_row["RT Window (min)"].values[0]
            rt, window = get_rt(rt_row, mode=mode)
            window = nstds * window
            window = max(min_window, min(window, max_window))
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time"] = rt
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time Window"] = window
            print(f"Updated RT for: {name} to {rt:.2f} +/- {window:.2f}")
        else:
            print(f"No RT found for: {name}")
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time Window"] = max_window
    return outdf



def tl_to_il(df):
    il_df = pd.DataFrame()
    il_df["Compound"] = df["Molecule Name"]
    il_df["m/z"] = df["Precursor Mz"]
    il_df["z"] = df["Precursor Charge"]

    rts = df["Explicit Retention Time"].to_numpy()
    rtwindows = df["Explicit Retention Time Window"].to_numpy()
    il_df["t start (min)"] = rts - (rtwindows / 2)
    b1 = il_df["t start (min)"] < 0
    il_df.loc[b1, "t start (min)"] = 0.1
    il_df["t stop (min)"] = rts + (rtwindows / 2)
    # il_df["Retention Time (min)"] = rts
    # il_df["RT Window (min)"] = rtwindows

    # il_df["HCD Collision Energy/Energies (%)"] = 30  # Default value, can be adjusted

    # Set polarity to be positive if charge is positive, negative if charge is negative
    il_df["Polarity"] = il_df["z"].apply(lambda x: "Positive" if x > 0 else "Negative")

    # Take abs of charge
    il_df["z"] = il_df["z"].abs()

    # Drop duplicates
    il_df = il_df.drop_duplicates().reset_index(drop=True)
    return il_df

def derplicate_il(df, mztol=0.2, rttol=1.0):
    good_rows = []
    df = df.sort_values(by=["m/z"]).reset_index(drop=True)
    df["Retention Time (min)"] = (df["t start (min)"] + df["t stop (min)"]) / 2
    # Group similar m/z and charge
    groups = []
    used = set()
    for i, row in df.iterrows():
        if i in used:
            continue
        group_indices = [i]
        for j in range(i + 1, len(df)):
            if j in used:
                continue
            other = df.loc[j]

            mz_close = abs(row["m/z"] - other["m/z"]) <= mztol

            same_charge = row["Polarity"] == other["Polarity"]

            rt_close = abs(
                row["Retention Time (min)"] -
                other["Retention Time (min)"]
            ) <= rttol

            if mz_close and same_charge and rt_close:
                group_indices.append(j)
        used.update(group_indices)
        group = df.loc[group_indices]
        groups.append(group)

    for group in groups:
        if len(group) == 1:
            good_rows.append(group.iloc[0].to_dict())
        else:
            # print(group.to_string())
            tstarts = group["t start (min)"].to_numpy()
            tstops = group["t stop (min)"].to_numpy()
            # Check if there are any overlaps
            overlap = False
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    if not (tstops[i] < tstarts[j] or tstops[j] < tstarts[i]):
                        overlap = True
                        break
                if overlap:
                    break
            if overlap:
                startt = min(tstarts)
                stopt = max(tstops)
                combined_row = {
                    "Compound": " & ".join(sorted(group["Compound"].astype(str).unique())),
                    "m/z": group["m/z"].mean(),
                    "z": group["z"].iloc[0],
                    "t start (min)": startt,
                    "t stop (min)": stopt,
                    # "HCD Collision Energy/Energies (%)": group["HCD Collision Energy/Energies (%)"].iloc[0]
                    "Polarity": group["Polarity"].iloc[0]
                }
                good_rows.append(combined_row)
                # print("Overlap found, new combined row:")
                # print(pd.Series(combined_row).to_string())
            else:
                for _, g in group.iterrows():
                    good_rows.append(g.to_dict())
                pass

    outdata = pd.DataFrame(good_rows).reset_index(drop=True)
    outdata["Retention Time (min)"] = (outdata["t start (min)"] + outdata["t stop (min)"]) / 2
    outdata["RT Window (min)"] = outdata["t stop (min)"] - outdata["t start (min)"]

    # Drop tstart and tstop columns
    outdata = outdata.drop(columns=["t start (min)", "t stop (min)"])

    return outdata

def simplify_class_name(name):
    if isinstance(name, str):
        name = name.replace("Ox", "")
        name = name.replace("Ether", "")
        if name.startswith("L"):
            name = name[1:]
        if name.startswith("ML") or name.startswith("DL"):
            name = name[2:]
        if name in hg_simplifier.keys():
            name = hg_simplifier[name]
    return name

def simplify_class_df(df, class_col="Ontology"):
    df[class_col] = df[class_col].apply(simplify_class_name)
    return df

def parse_spec_string(spec_string, threshold=0, norm=False):
    # Parse a string of the form "m/z:intensity m/z:intensity,..." into a numpy array
    if pd.isna(spec_string):
        return np.array([]).reshape(0, 2)
    pairs = spec_string.split(" ")
    data = []
    for pair in pairs:
        try:
            mz, intensity = pair.split(":")
            data.append([float(mz), float(intensity)])
        except ValueError:
            continue
    data = np.array(data)
    # Apply threshold
    data = data[data[:, 1] > threshold * np.amax(data[:, 1])]

    if norm:
        data[:, 1] = data[:, 1] / np.amax(data[:, 1])
    return data

def match_frag_to_spec(spec_data, frag_mz, mz_tol=0.05):
    # spec_data is a numpy array of shape (n, 2) with columns [m/z, intensity]
    # frag_mz is a float
    # mz_tol is a float
    if spec_data.size == 0:
        return None
    mz_diff = np.abs(spec_data[:, 0] - frag_mz)
    min_diff_index = np.argmin(mz_diff)
    if mz_diff[min_diff_index] <= mz_tol:
        return True
    else:
        return False

def spec_to_nl(spec_data, precursor_mz):
    # spec_data is a numpy array of shape (n, 2) with columns [m/z, intensity]
    # precursor_mz is a float
    if spec_data.size == 0:
        return np.array([]).reshape(0, 2)
    nl_data = np.copy(spec_data)
    nl_data[:, 0] = precursor_mz - spec_data[:, 0]
    return nl_data

def mzrt_plot(df, class_col="Ontology", title="", simplify=True, write_file=False, labeln=False, rtcol="Average Rt(min)", mzcol="Average Mz"):
    # ------------------------------
    # Plot setup
    # ------------------------------
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 18

    if simplify:
        df = simplify_class_df(df, class_col=class_col)
    else:
        # Apply hg_simplifier to class_col for better color mapping
        df[class_col] = df[class_col].apply(lambda x: hg_simplifier.get(x) if x in hg_simplifier.keys() else x)

    # Legend order: alphabetical, but put Other last if present
    unique_classes = sorted([c for c in df[class_col].unique() if c != "Other"])
    if "Other" in df[class_col].unique():
        unique_classes.append("Other")

    # Build legend
    custom_legend = [
         mpl.lines.Line2D([0], [0], marker='o', color=get_color_for_class(cls), label=cls,
               markerfacecolor=get_color_for_class(cls), linestyle='None', markersize=10)
        for cls in unique_classes
    ]

    # Set Color column for plotting
    df["Color"] = df[class_col].apply(get_color_for_class)


    # ------------------------------
    # Plot
    # ------------------------------
    plt.figure(figsize=(12, 8))
    plt.scatter(df[rtcol], df[mzcol], c=df["Color"], s=80)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel("Retention Time (min)")
    plt.ylabel("m/z", fontstyle="italic")
    plt.title(title)
    plt.xticks(np.arange(0.0, 18, 2.5))
    plt.yticks(np.arange(0.0, 1600.1, 200))

    if labeln:
        # Add total number of lipids on the plot
        total_lipids = len(df)
        plt.text(
            0.98, 0.98, f"Total lipids: {total_lipids}",
            horizontalalignment='right',
            verticalalignment='top',
            transform=plt.gca().transAxes,
            fontsize=20,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=5)
        )

    ncol = min(6, len(unique_classes))
    plt.legend(handles=custom_legend, loc='upper center', bbox_to_anchor=(0.5, -0.1),
               ncol=ncol, frameon=False)

    plt.subplots_adjust(bottom=0.25)
    plt.tight_layout()
    if write_file:
        plt.savefig(f"{title}_mzrt_plot.pdf", bbox_inches='tight')
        plt.savefig(f"{title}_mzrt_plot.png", bbox_inches='tight', dpi=300)

    plt.show()

def write_colored_excel(df, column, filename, color_map=None):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
        # Convert the dataframe to an XlsxWriter Excel object.
        df.to_excel(writer, sheet_name='Sheet1', index=False)

        # Get the xlsxwriter workbook and worksheet objects.
        workbook = writer.book
        worksheet = writer.sheets['Sheet1']

        # Define a default color map if none is provided
        if color_map is None:
            unique_values = df[column].unique()
            colors = plt.cm.get_cmap('tab20', len(unique_values)).colors
            color_map = {val: colors[i] for i, val in enumerate(unique_values)}

        # Apply formatting based on the column values
        for row_num, value in enumerate(df[column], start=1):  # start=1 to skip header
            color = color_map.get(value, (1, 1, 1))  # Default to white if not found
            color = mpl.colors.to_hex(color)
            hex_color = '#%02x%02x%02x' % (int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16))
            cell_format = workbook.add_format({'bg_color': hex_color})
            worksheet.write(row_num, df.columns.get_loc(column), value, cell_format)


def get_tail_carbons(tail):
    """
    Extracts carbon count from tail strings like:
    18:1;O2
    16:0
    24:1;O
    """
    if pd.isna(tail):
        return np.nan

    tail = str(tail).strip()

    if tail == "":
        return np.nan

    match = re.match(r"^(\d+):", tail)

    if match is None:
        return np.nan

    return float(match.group(1))


def get_tail_unsaturation(tail):
    """
    Extracts unsaturation count from tail strings like:
    18:1;O2
    16:0
    24:1;O
    """
    if pd.isna(tail):
        return np.nan

    tail = str(tail).strip()

    if tail == "":
        return np.nan

    match = re.match(r"^\d+:(\d+)", tail)

    if match is None:
        return np.nan

    return float(match.group(1))


def flag_weird_lipids(
    df,
    lipid_col="Metabolite name",
    class_col="Ontology",
    adduct_col="Adduct type",

    # General total-unsaturation limits by class
    max_unsaturation_by_class=None,

    # SM-specific structural limits
    sm_min_backbone_carbons=16,
    sm_max_backbone_carbons=20,
    sm_max_acyl_carbons=26,
    sm_max_backbone_unsaturation=2,
    sm_max_acyl_unsaturation=3,
    sm_min_total_carbons=30
):
    """
    Applies general and class-specific lipid sanity checks.

    General unsaturation rule:
        Classes listed in max_unsaturation_by_class are removed when their
        total number of double bonds exceeds the class-specific limit.

    SM-specific rules:
        Speciated SMs:
            - T1 backbone carbon limits
            - T2 acyl-tail carbon limit
            - T1 backbone unsaturation limit
            - T2 acyl-tail unsaturation limit

        Non-speciated SMs:
            - Minimum total carbon count

    PC-specific rule:
        - Flags PCs found only in positive mode and not in negative mode.

    No helper or flag columns are added to the returned DataFrame.
    """

    df = df.copy()

    if max_unsaturation_by_class is None:
        max_unsaturation_by_class = {
            "PC": 6,
            "SM": 4
        }

    if "Overall Match Status" not in df.columns:
        df["Overall Match Status"] = True

    if "Matched Required Fragments" not in df.columns:
        df["Matched Required Fragments"] = ""

    df["Matched Required Fragments"] = (
        df["Matched Required Fragments"]
        .fillna("")
        .astype(str)
    )

    clean_class = (
        df[class_col]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    clean_lipid_name = (
        df[lipid_col]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    cleaned_adducts = df[adduct_col].apply(
        lambda value: fix_adduct(value, 1, 0)
    )

    identification_polarity = pd.Series(
        np.where(
            cleaned_adducts.astype(str).str.endswith("+"),
            "Positive",
            np.where(
                cleaned_adducts.astype(str).str.endswith("-"),
                "Negative",
                "Unknown"
            )
        ),
        index=df.index
    )

    if "T1" in df.columns:
        t1_carbons = df["T1"].apply(get_tail_carbons)
        t1_unsaturation = df["T1"].apply(get_tail_unsaturation)
    else:
        t1_carbons = pd.Series(np.nan, index=df.index)
        t1_unsaturation = pd.Series(np.nan, index=df.index)

    if "T2" in df.columns:
        t2_carbons = df["T2"].apply(get_tail_carbons)
        t2_unsaturation = df["T2"].apply(get_tail_unsaturation)
    else:
        t2_carbons = pd.Series(np.nan, index=df.index)
        t2_unsaturation = pd.Series(np.nan, index=df.index)

    # For speciated lipids, add T1 and T2 unsaturation.
    # For sum-composition lipids, T2 is missing and T1 contains the total.
    total_unsaturation = pd.Series(
        np.where(
            t2_unsaturation.notna(),
            t1_unsaturation.fillna(0)
            + t2_unsaturation.fillna(0),
            t1_unsaturation
        ),
        index=df.index
    )

    # ==============================================================
    # GENERAL TOTAL-UNSATURATION CHECK
    # ==============================================================

    excessive_unsaturation = pd.Series(False, index=df.index)

    print("\nGeneral lipid unsaturation checks:")

    for lipid_class, max_unsaturation in max_unsaturation_by_class.items():
        class_excessive_unsaturation = (
            clean_class.eq(lipid_class)
            & total_unsaturation.notna()
            & (total_unsaturation > max_unsaturation)
        )

        excessive_unsaturation = (
            excessive_unsaturation
            | class_excessive_unsaturation
        )

        df.loc[
            class_excessive_unsaturation,
            "Matched Required Fragments"
        ] += (
            f" FLAGGED: {lipid_class} total double bonds "
            f"> {max_unsaturation};"
        )

        print(
            f"{lipid_class} rows removed for more than "
            f"{max_unsaturation} double bonds:",
            int(class_excessive_unsaturation.sum())
        )

    # ==============================================================
    # SM-SPECIFIC CHECKS
    # ==============================================================

    sm_mask = clean_class.eq("SM")

    speciated_sm_mask = (
        sm_mask
        & clean_lipid_name.str.contains(r"\|", na=False)
    )

    sum_sm_mask = (
        sm_mask
        & ~clean_lipid_name.str.contains(r"\|", na=False)
    )

    sm_small_backbone = (
        speciated_sm_mask
        & t1_carbons.notna()
        & (t1_carbons < sm_min_backbone_carbons)
    )

    sm_large_backbone = (
        speciated_sm_mask
        & t1_carbons.notna()
        & (t1_carbons > sm_max_backbone_carbons)
    )

    sm_large_acyl_tail = (
        speciated_sm_mask
        & t2_carbons.notna()
        & (t2_carbons > sm_max_acyl_carbons)
    )

    sm_unsaturated_backbone = (
        speciated_sm_mask
        & t1_unsaturation.notna()
        & (t1_unsaturation > sm_max_backbone_unsaturation)
    )

    sm_unsaturated_acyl_tail = (
        speciated_sm_mask
        & t2_unsaturation.notna()
        & (t2_unsaturation > sm_max_acyl_unsaturation)
    )

    sm_small_sum_composition = (
        sum_sm_mask
        & t1_carbons.notna()
        & (t1_carbons < sm_min_total_carbons)
    )

    sm_weird_lipid = (
        sm_small_backbone
        | sm_large_backbone
        | sm_large_acyl_tail
        | sm_unsaturated_backbone
        | sm_unsaturated_acyl_tail
        | sm_small_sum_composition
    )

    df.loc[
        sm_small_backbone,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: speciated SM backbone T1 carbon count "
        f"< {sm_min_backbone_carbons};"
    )

    df.loc[
        sm_large_backbone,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: speciated SM backbone T1 carbon count "
        f"> {sm_max_backbone_carbons};"
    )

    df.loc[
        sm_large_acyl_tail,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: speciated SM acyl tail T2 carbon count "
        f"> {sm_max_acyl_carbons};"
    )

    df.loc[
        sm_unsaturated_backbone,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: speciated SM backbone T1 unsaturation "
        f"> {sm_max_backbone_unsaturation};"
    )

    df.loc[
        sm_unsaturated_acyl_tail,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: speciated SM acyl tail T2 unsaturation "
        f"> {sm_max_acyl_unsaturation};"
    )

    df.loc[
        sm_small_sum_composition,
        "Matched Required Fragments"
    ] += (
        f" FLAGGED: non-speciated SM total carbon count "
        f"< {sm_min_total_carbons};"
    )

    # ==============================================================
    # PC-SPECIFIC POLARITY CHECK
    # ==============================================================

    pc_mask = clean_class.eq("PC")

    pc_polarity_by_name = (
        pd.DataFrame({
            "_Clean_Lipid_Name": clean_lipid_name.loc[pc_mask],
            "_Polarity": identification_polarity.loc[pc_mask]
        })
        .groupby("_Clean_Lipid_Name")["_Polarity"]
        .agg(lambda values: set(values))
    )

    positive_only_pc_names = pc_polarity_by_name[
        pc_polarity_by_name.apply(
            lambda polarities:
                "Positive" in polarities
                and "Negative" not in polarities
        )
    ].index

    pc_positive_only = (
        pc_mask
        & clean_lipid_name.isin(positive_only_pc_names)
    )

    df.loc[
        pc_positive_only,
        "Matched Required Fragments"
    ] += (
        " FLAGGED: PC identified only in positive mode; "
        "no negative-mode identification found;"
    )

    weird_lipid = (
        excessive_unsaturation
        | sm_weird_lipid
        | pc_positive_only
    )

    df.loc[
        weird_lipid,
        "Overall Match Status"
    ] = False

    print("\nClass-specific weird lipid checks:")

    print(
        "SM rows flagged by structural rules:",
        int(sm_weird_lipid.sum())
    )

    print(
        "PC rows flagged as positive-mode only:",
        int(pc_positive_only.sum())
    )

    print(
        "Unique positive-mode-only PC names:",
        len(positive_only_pc_names)
    )

    print(
        "Total rows flagged by all lipid sanity checks:",
        int(weird_lipid.sum())
    )

    # Remove all listed classes that exceed their unsaturation limit.
    df = df[~excessive_unsaturation].copy()

    return df

def plot_rejected_quality_density(
    candidate_df,
    output_basename="Rejected_Quality_Density"
):
    plot_df = candidate_df.copy()

    accepted_values = [
        True, "TRUE", "True", "true", 1, "1"
    ]

    # Keep rejected candidates only
    plot_df = plot_df[
        ~plot_df["Overall Match Status"].isin(accepted_values)
    ].copy()

    fragment_text = (
        plot_df["Matched Required Fragments"]
        .fillna("")
        .astype(str)
    )

    plot_df["Fragments Found"] = (
        fragment_text.str.count(r":\s*Yes")
    )

    plot_df["Fragments Missing"] = (
        fragment_text.str.count(r":\s*No")
    )

    plot_df["Fragments Evaluated"] = (
        plot_df["Fragments Found"]
        + plot_df["Fragments Missing"]
    )

    plot_df["Fragment Match Fraction"] = np.where(
        plot_df["Fragments Evaluated"] > 0,
        (
            plot_df["Fragments Found"]
            / plot_df["Fragments Evaluated"]
        ),
        np.nan
    )

    plot_df["Reverse dot product"] = pd.to_numeric(
        plot_df["Reverse dot product"],
        errors="coerce"
    )

    plot_df = plot_df.dropna(
        subset=[
            "Reverse dot product",
            "Fragment Match Fraction"
        ]
    )

    fig, ax = plt.subplots(figsize=(9, 7))

    density = ax.hist2d(
        plot_df["Reverse dot product"],
        plot_df["Fragment Match Fraction"],
        bins=[20, 10],
        cmap="viridis"
    )

    fig.colorbar(
        density[3],
        ax=ax,
        label="Number of Rejected Candidates"
    )

    ax.set_xlabel("Reverse Dot Product")
    ax.set_ylabel("Fraction of Required Fragments Detected")

    ax.set_title(
        "Common Raw-Data Characteristics of Rejected Lipid Candidates"
    )

    ax.set_ylim(0, 1)

    fig.tight_layout()

    fig.savefig(
        f"{output_basename}.png",
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        f"{output_basename}.pdf",
        bbox_inches="tight"
    )

    return plot_df

def plot_class_dataset_heatmap(
    df,
    datasets,
    output_basename="Class_Dataset_Heatmap"
):
    """
    Plot the number of unique confirmed lipid IDs in each lipid class
    for each dataset.
    """

    plot_df = df.copy()

    # Clean dataset and lipid-class names
    plot_df["Dataset"] = (
        plot_df["Dataset"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df["Ontology"] = (
        plot_df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    # Keep confirmed IDs only
    plot_df = plot_df[
        plot_df["Overall Match Status"].isin(
            [True, "TRUE", "True", "true", 1, "1"]
        )
    ].copy()

    # Keep only the selected datasets
    plot_df = plot_df[
        plot_df["Dataset"].isin(datasets)
    ].copy()

    # Match speciated and unspeciated versions using sum composition
    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    # Count each lipid only once per class and dataset
    plot_df = plot_df.drop_duplicates(
        subset=[
            "Dataset",
            "Ontology",
            "Lipid Match Key"
        ]
    )

    # Create class-by-dataset count matrix
    heatmap_df = (
        plot_df
        .groupby(
            ["Ontology", "Dataset"]
        )
        .size()
        .unstack(fill_value=0)
        .reindex(
            columns=datasets,
            fill_value=0
        )
    )

    # Sort classes by total number of IDs
    heatmap_df["Total"] = heatmap_df.sum(axis=1)

    heatmap_df = (
        heatmap_df
        .sort_values(
            "Total",
            ascending=False
        )
        .drop(columns=["Total"])
    )

    print("\nConfirmed lipid counts by class and dataset:")
    print(heatmap_df.to_string())

    # Create heatmap
    fig, ax = plt.subplots(
        figsize=(
            8,
            max(6, len(heatmap_df) * 0.4)
        )
    )

    image = ax.imshow(
        heatmap_df.values,
        aspect="auto",
        cmap="viridis"
    )

    ax.set_xticks(
        range(len(heatmap_df.columns))
    )

    ax.set_xticklabels(
        heatmap_df.columns,
        rotation=45,
        ha="right"
    )

    ax.set_yticks(
        range(len(heatmap_df.index))
    )

    ax.set_yticklabels(
        heatmap_df.index
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Lipid Class")

    ax.set_title(
        "Unique Confirmed Lipid IDs by Class and Dataset"
    )

    # Print the raw count inside every cell
    for row_index in range(len(heatmap_df.index)):
        for col_index in range(len(heatmap_df.columns)):

            value = int(
                heatmap_df.iloc[
                    row_index,
                    col_index
                ]
            )

            ax.text(
                col_index,
                row_index,
                str(value),
                ha="center",
                va="center",
                fontsize=8
            )

    fig.colorbar(
        image,
        ax=ax,
        label="Unique Confirmed Lipid IDs"
    )

    fig.tight_layout()

    png_file = f"{output_basename}.png"
    pdf_file = f"{output_basename}.pdf"

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    print("\nHeatmap saved as:")
    print(os.path.abspath(png_file))
    print(os.path.abspath(pdf_file))

    return heatmap_df

def plot_class_consistency(
    df,
    datasets,
    output_basename="Class_Consistency"
):
    """
    For each lipid class, calculate the average percentage of datasets
    in which its confirmed unique lipid identities were detected.
    """

    plot_df = df.copy()

    # Keep confirmed identifications only
    plot_df = plot_df[
        plot_df["Overall Match Status"].isin(
            [True, "TRUE", "True", 1, "1"]
        )
    ].copy()

    # Clean relevant columns
    plot_df["Dataset"] = (
        plot_df["Dataset"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df["Ontology"] = (
        plot_df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    # Convert detailed names to sum-composition names
    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    # Keep requested datasets
    plot_df = plot_df[
        plot_df["Dataset"].isin(datasets)
    ].copy()

    # One occurrence of each lipid per dataset
    plot_df = plot_df.drop_duplicates(
        subset=[
            "Dataset",
            "Ontology",
            "Lipid Match Key"
        ]
    )

    # Number of datasets containing each lipid
    lipid_coverage = (
        plot_df
        .groupby(
            ["Ontology", "Lipid Match Key"]
        )["Dataset"]
        .nunique()
        .reset_index(name="Datasets Present")
    )

    lipid_coverage["Consistency Percent"] = (
        lipid_coverage["Datasets Present"]
        / len(datasets)
        * 100
    )

    # Average consistency for each class
    class_summary = (
        lipid_coverage
        .groupby("Ontology")
        .agg(
            Mean_Consistency=(
                "Consistency Percent",
                "mean"
            ),
            Number_of_Lipids=(
                "Lipid Match Key",
                "nunique"
            )
        )
        .reset_index()
        .sort_values(
            "Mean_Consistency",
            ascending=False
        )
    )

    print("\nClass consistency values:")
    print(class_summary.to_string(index=False))

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.bar(
        class_summary["Ontology"],
        class_summary["Mean_Consistency"],
        edgecolor="black"
    )

    ax.set_xlabel("Lipid Class")
    ax.set_ylabel("Mean Dataset Consistency (%)")
    ax.set_title("Consistency of Lipid Classes Across Datasets")
    ax.set_ylim(0, 100)

    for i, value in enumerate(
        class_summary["Mean_Consistency"]
    ):
        ax.text(
            i,
            value,
            f"{value:.1f}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    png_file = f"{output_basename}.png"
    pdf_file = f"{output_basename}.pdf"

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    print("\nClass-consistency plots saved as:")
    print(os.path.abspath(png_file))
    print(os.path.abspath(pdf_file))

    plt.show()
    plt.close(fig)

    return class_summary


def plot_total_ids_by_dataset(
    df,
    datasets,
    dataset_col="Dataset",
    bar_color=None,
    output_basename="Confirmed_Lipid_IDs_by_Dataset"
):
    """
    Plot one confirmed sum-composition lipid ID per class and dataset.

    Speciated and unspeciated names with the same sum composition are
    counted as one identity within a dataset.
    """
    plot_df = df.copy()

    plot_df[dataset_col] = (
        plot_df[dataset_col]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df["Ontology"] = (
        plot_df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    # Keep confirmed identifications only.
    plot_df = plot_df[
        plot_df["Overall Match Status"].isin(
            [True, "TRUE", "True", 1, "1"]
        )
    ].copy()

    # Convert detailed names to sum-composition names.
    plot_df["Plot Lipid Name"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    # Keep requested datasets and count each lipid once per dataset.
    plot_df = plot_df[
        plot_df[dataset_col].isin(datasets)
    ].drop_duplicates(
        subset=[
            dataset_col,
            "Ontology",
            "Plot Lipid Name"
        ]
    )

    counts = (
        plot_df
        .groupby(dataset_col)
        .size()
        .reindex(datasets, fill_value=0)
    )

    print("\nCounts used for confirmed-ID plot:")
    print(counts.to_string())

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(
        counts.index,
        counts.values,
        color=bar_color,
        edgecolor="black"
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Unique Confirmed Lipid IDs")
    ax.set_title("Confirmed Lipid IDs by Dataset")

    for i, value in enumerate(counts.values):
        ax.text(
            i,
            value,
            str(int(value)),
            ha="center",
            va="bottom"
        )

    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    png_file = f"{output_basename}.png"
    pdf_file = f"{output_basename}.pdf"

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    print("\nPlots saved as:")
    print(os.path.abspath(png_file))
    print(os.path.abspath(pdf_file))

    plt.show(block=True)
    plt.close(fig)


def export_lipid_dataset_presence(
    df,
    datasets,
    lipid_classes,
    lipid_col="Lipid",
    class_col="Ontology",
    dataset_col="Dataset",
    output_file="Lipid_Dataset_Presence.xlsx"
):
    """
    Exports dataset presence/absence information for specified lipid classes.

    Output includes:
        - One row per unique lipid
        - 1/0 presence columns for each dataset
        - Number of datasets containing the lipid
        - Number of datasets missing the lipid
        - Names of datasets containing the lipid
        - Names of datasets missing the lipid
        - Separate Excel sheets grouped by number of datasets present
    """

    presence_df = df.copy()

    presence_df = presence_df[
        presence_df["Overall Match Status"].isin(
            [True, "TRUE", "True", 1, "1"]
        )
    ].copy()

    # Clean relevant text columns
    for col in [lipid_col, class_col, dataset_col]:
        presence_df[col] = (
            presence_df[col]
            .astype(str)
            .str.replace("\u00a0", " ", regex=False)
            .str.strip()
        )
    presence_df["Lipid Match Key"] = (
        presence_df[lipid_col]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    presence_df["Is Speciated"] = (
        presence_df[lipid_col]
        .astype(str)
        .str.contains("|", regex=False)
    )
    # Keep only requested datasets and lipid classes
    presence_df = presence_df[
        presence_df[dataset_col].isin(datasets)
        & presence_df[class_col].isin(lipid_classes)
    ].copy()

    # Remove rows without a meaningful lipid name
    presence_df = presence_df[
        presence_df[lipid_col].notna()
        & ~presence_df[lipid_col].isin(["", "nan", "None"])
    ]

    # One row per lipid/class/dataset before making the table
    presence_df = presence_df[
        [
            class_col,
            lipid_col,
            "Lipid Match Key",
            dataset_col,
            "Is Speciated"
        ]
    ].drop_duplicates()

    # Build presence/absence matrix
    presence_matrix = (
        presence_df
        .assign(Present=1)
        .pivot_table(
            index=[class_col, "Lipid Match Key"],
            columns=dataset_col,
            values="Present",
            aggfunc="max",
            fill_value=0
        )
        .reindex(columns=datasets, fill_value=0)
        .reset_index()
    )

    # Prefer a speciated name for display
    display_names = (
        presence_df
        .sort_values("Is Speciated", ascending=False)
        .drop_duplicates([class_col, "Lipid Match Key"])
        [[class_col, "Lipid Match Key", lipid_col]]
    )

    presence_matrix = presence_matrix.merge(
        display_names,
        on=[class_col, "Lipid Match Key"],
        how="left"
    )

    # List datasets where only the unspeciated form was found
    speciation_status = (
        presence_df
        .groupby(
            [class_col, "Lipid Match Key", dataset_col]
        )["Is Speciated"]
        .any()
        .reset_index()
    )

    not_speciated = (
        speciation_status[
            speciation_status["Is Speciated"] == False
            ]
        .groupby(
            [class_col, "Lipid Match Key"]
        )[dataset_col]
        .apply(lambda x: ", ".join(sorted(x)))
        .reset_index(name="Not Speciated In")
    )

    presence_matrix = presence_matrix.merge(
        not_speciated,
        on=[class_col, "Lipid Match Key"],
        how="left"
    )

    presence_matrix["Not Speciated In"] = (
        presence_matrix["Not Speciated In"].fillna("")
    )

    # Count datasets present and missing
    presence_matrix["Datasets Present"] = (
        presence_matrix[datasets].sum(axis=1)
    )

    presence_matrix["Datasets Missing"] = (
        len(datasets) - presence_matrix["Datasets Present"]
    )

    # List dataset names
    presence_matrix["Present In"] = presence_matrix.apply(
        lambda row: ", ".join(
            dataset for dataset in datasets if row[dataset] == 1
        ),
        axis=1
    )

    presence_matrix["Missing From"] = presence_matrix.apply(
        lambda row: ", ".join(
            dataset for dataset in datasets if row[dataset] == 0
        ),
        axis=1
    )

    # Sort by class, then by greatest coverage, then lipid name
    presence_matrix = presence_matrix.sort_values(
        by=[class_col, "Datasets Present", lipid_col],
        ascending=[True, False, True]
    )

    # Summary by class and number of datasets present
    summary = (
        presence_matrix
        .groupby([class_col, "Datasets Present"])
        .size()
        .reset_index(name="Number of Lipids")
        .sort_values([class_col, "Datasets Present"])
    )

    presence_matrix = presence_matrix.drop(
        columns=["Lipid Match Key"]
    )

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:

        # Complete table
        presence_matrix.to_excel(
            writer,
            sheet_name="All Lipids",
            index=False
        )

        # Summary counts
        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False
        )

        # Separate sheet for each presence count
        for number_present in range(len(datasets) + 1):
            subset = presence_matrix[
                presence_matrix["Datasets Present"] == number_present
            ]

            if not subset.empty:
                subset.to_excel(
                    writer,
                    sheet_name=f"Present in {number_present}",
                    index=False
                )

        # Separate sheets for each requested class
        for lipid_class in lipid_classes:
            class_subset = presence_matrix[
                presence_matrix[class_col] == lipid_class
            ]

            if not class_subset.empty:
                safe_sheet_name = str(lipid_class)[:31]

                class_subset.to_excel(
                    writer,
                    sheet_name=safe_sheet_name,
                    index=False
                )

    print(f"\nPresence/absence spreadsheet saved as:")
    print(os.path.abspath(output_file))

    print("\nNumber of lipids by class and dataset coverage:")
    print(summary.to_string(index=False))

    return presence_matrix

def plot_removed_ids_by_class(
    removed_df,
    output_basename="Removed_IDs_by_Class"
):
    plot_df = removed_df.copy()

    plot_df["Ontology"] = (
        plot_df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    # Count each removed lipid only once across all datasets
    plot_df = plot_df.drop_duplicates(
        subset=[
            "Ontology",
            "Lipid Match Key"
        ]
    )

    removed_counts = (
        plot_df["Ontology"]
        .value_counts()
        .sort_values(ascending=False)
    )

def plot_class_coverage_distribution(
    df,
    datasets,
    output_basename="Class_Coverage_Distribution"
):
    plot_df = df.copy()

    plot_df["Dataset"] = (
        plot_df["Dataset"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df["Ontology"] = (
        plot_df["Ontology"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    plot_df = plot_df[
        plot_df["Overall Match Status"].isin(
            [True, "TRUE", "True", 1, "1"]
        )
    ].copy()

    plot_df = plot_df[
        plot_df["Dataset"].isin(datasets)
    ].copy()

    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    plot_df = plot_df.drop_duplicates(
        subset=[
            "Dataset",
            "Ontology",
            "Lipid Match Key"
        ]
    )

    lipid_coverage = (
        plot_df
        .groupby(
            ["Ontology", "Lipid Match Key"]
        )["Dataset"]
        .nunique()
        .reset_index(name="Datasets Present")
    )

    coverage_counts = (
        lipid_coverage
        .groupby(
            ["Ontology", "Datasets Present"]
        )
        .size()
        .unstack(fill_value=0)
        .reindex(columns=range(1, len(datasets) + 1), fill_value=0)
    )

    coverage_percent = (
        coverage_counts
        .div(coverage_counts.sum(axis=1), axis=0)
        * 100
    )

    # Sort classes by percentage found in 3 or more datasets
    coverage_percent["3+ Total"] = (
        coverage_percent[
            [col for col in coverage_percent.columns if col >= 3]
        ]
        .sum(axis=1)
    )

    coverage_percent = coverage_percent.sort_values(
        "3+ Total",
        ascending=True
    )

    coverage_percent = coverage_percent.drop(
        columns=["3+ Total"]
    )

    print("\nCoverage distribution by class (%):")
    print(coverage_percent.to_string())

    fig, ax = plt.subplots(
        figsize=(10, max(6, len(coverage_percent) * 0.35))
    )

    cmap = plt.get_cmap("viridis")

    coverage_colors = [
        cmap(i / (len(datasets) - 1))
        for i in range(len(datasets))
    ]

    left = np.zeros(len(coverage_percent))

    for coverage_number in range(1, len(datasets) + 1):
        values = coverage_percent[coverage_number].values

        ax.barh(
            coverage_percent.index,
            values,
            left=left,
            label=(
                f"{coverage_number} dataset"
                if coverage_number == 1
                else f"{coverage_number} datasets"
            ),
            edgecolor="black",
            color=coverage_colors[coverage_number - 1]
        )

        left += values
    ax.set_xlabel("Unique Lipids in Each Coverage Group (%)")
    ax.set_ylabel("Lipid Class")
    ax.set_title("Dataset Coverage Distribution by Lipid Class")
    ax.set_xlim(0, 100)

    ax.legend(
        title="Detected In",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    fig.tight_layout()

    fig.savefig(
        f"{output_basename}.png",
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        f"{output_basename}.pdf",
        bbox_inches="tight"
    )

    return coverage_percent

def plot_true_false_quality_distributions(
    candidate_df,
    output_basename="True_False_Quality_Distributions"
):
    plot_df = candidate_df.copy()

    true_values = [True, "TRUE", "True", "true", 1, "1"]

    plot_df["Confirmed"] = (
        plot_df["Overall Match Status"]
        .isin(true_values)
    )

    fragment_text = (
        plot_df["Matched Required Fragments"]
        .fillna("")
        .astype(str)
    )

    plot_df["Fragments Found"] = (
        fragment_text.str.count(r":\s*Yes")
    )

    plot_df["Fragments Missing"] = (
        fragment_text.str.count(r":\s*No")
    )

    plot_df["Fragments Evaluated"] = (
        plot_df["Fragments Found"]
        + plot_df["Fragments Missing"]
    )

    plot_df["Fragment Match Fraction"] = np.where(
        plot_df["Fragments Evaluated"] > 0,
        (
            plot_df["Fragments Found"]
            / plot_df["Fragments Evaluated"]
        ),
        np.nan
    )

    plot_df["Reverse dot product"] = pd.to_numeric(
        plot_df["Reverse dot product"],
        errors="coerce"
    )

    features = [
        (
            "Reverse dot product",
            "Reverse Dot Product",
            np.linspace(0, 1, 25)
        ),
        (
            "Fragment Match Fraction",
            "Fraction of Required Fragments Matched",
            np.linspace(0, 1, 12)
        )
    ]

    fig, axes = plt.subplots(
        ncols=2,
        figsize=(10, 5)
    )

    for ax, (column, title, bins) in zip(
        axes,
        features
    ):
        false_values = (
            plot_df.loc[
                ~plot_df["Confirmed"],
                column
            ]
            .dropna()
        )

        true_values_data = (
            plot_df.loc[
                plot_df["Confirmed"],
                column
            ]
            .dropna()
        )

        ax.hist(
            false_values,
            bins=bins,
            alpha=0.65,
            label="False",
            density=True
        )

        ax.hist(
            true_values_data,
            bins=bins,
            alpha=0.65,
            label="True",
            density=True
        )

        ax.set_title(title)
        ax.set_xlabel(title)
        ax.set_ylabel("Density")
        ax.legend()

    fig.suptitle(
        "Raw-Data Characteristics of Confirmed and Rejected Lipid IDs"
    )

    fig.tight_layout()

    fig.savefig(
        f"{output_basename}.png",
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        f"{output_basename}.pdf",
        bbox_inches="tight"
    )

    return plot_df

def plot_quality_feature_separation(
    candidate_df,
    output_basename="Quality_Feature_Separation"
):
    plot_df = candidate_df.copy()

    true_values = [True, "TRUE", "True", "true", 1, "1"]

    plot_df["Confirmed"] = (
        plot_df["Overall Match Status"]
        .isin(true_values)
    )

    fragment_text = (
        plot_df["Matched Required Fragments"]
        .fillna("")
        .astype(str)
    )

    plot_df["Fragments Found"] = (
        fragment_text.str.count(r":\s*Yes")
    )

    plot_df["Fragments Missing"] = (
        fragment_text.str.count(r":\s*No")
    )

    plot_df["Fragments Evaluated"] = (
        plot_df["Fragments Found"]
        + plot_df["Fragments Missing"]
    )

    plot_df["Fragment Match Fraction"] = np.where(
        plot_df["Fragments Evaluated"] > 0,
        (
            plot_df["Fragments Found"]
            / plot_df["Fragments Evaluated"]
        ),
        np.nan
    )

    feature_columns = [
        "Reverse dot product",
        "Fragment Match Fraction"
    ]
    feature_labels = {
        "Reverse dot product": "Reverse dot product",
        "Fragment Match Fraction": "Fragment-match fraction"
    }
    separation_scores = {}

    for feature in feature_columns:
        plot_df[feature] = pd.to_numeric(
            plot_df[feature],
            errors="coerce"
        )

        true_group = (
            plot_df.loc[
                plot_df["Confirmed"],
                feature
            ]
            .dropna()
        )

        false_group = (
            plot_df.loc[
                ~plot_df["Confirmed"],
                feature
            ]
            .dropna()
        )

        pooled_sd = np.sqrt(
            (
                true_group.var(ddof=1)
                + false_group.var(ddof=1)
            )
            / 2
        )

        if pooled_sd == 0 or np.isnan(pooled_sd):
            score = 0
        else:
            score = abs(
                true_group.mean()
                - false_group.mean()
            ) / pooled_sd

        separation_scores[
            feature_labels[feature]
        ] = score

    score_series = (
        pd.Series(separation_scores)
        .sort_values(ascending=False)
    )

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.bar(
        score_series.index,
        score_series.values,
        edgecolor="black"
    )

    ax.set_xlabel("Raw-Data Feature")
    ax.set_ylabel("Standardized Difference")
    ax.set_title(
        "Features That Best Separate Confirmed and Rejected IDs"
    )

    ax.tick_params(
        axis="x",
        rotation=45
    )

    fig.tight_layout()

    fig.savefig(
        f"{output_basename}.png",
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        f"{output_basename}.pdf",
        bbox_inches="tight"
    )

    return score_series

def plot_rejected_ids_by_dataset(
    candidate_df,
    datasets,
    bar_color=None,
    output_basename="Rejected_Lipid_IDs_by_Dataset"
):
    """
    Plot the number of unique lipid identities rejected by fragment
    validation in each dataset.

    Detailed and unspeciated names sharing the same sum composition
    are counted as one lipid within a dataset.
    """

    plot_df = candidate_df.copy()

    # Clean dataset names
    plot_df["Dataset"] = (
        plot_df["Dataset"]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    # Convert Overall Match Status into a reliable Boolean
    confirmed_values = [
        True,
        "TRUE",
        "True",
        "true",
        1,
        "1"
    ]

    plot_df["Confirmed"] = (
        plot_df["Overall Match Status"]
        .isin(confirmed_values)
    )

    # Keep rejected IDs only
    plot_df = plot_df[
        ~plot_df["Confirmed"]
    ].copy()

    # Keep the selected datasets
    plot_df = plot_df[
        plot_df["Dataset"].isin(datasets)
    ].copy()

    # Convert detailed names to sum-composition identities
    plot_df["Lipid Match Key"] = (
        plot_df["Metabolite name"]
        .astype(str)
        .str.split("|", n=1)
        .str[0]
        .str.strip()
    )

    # Remove invalid names
    plot_df = plot_df[
        ~plot_df["Lipid Match Key"].isin(
            ["", "nan", "None"]
        )
    ].copy()

    # Count each rejected lipid once per dataset
    plot_df = plot_df.drop_duplicates(
        subset=[
            "Dataset",
            "Ontology",
            "Lipid Match Key"
        ]
    )

    rejected_counts = (
        plot_df
        .groupby("Dataset")
        .size()
        .reindex(datasets, fill_value=0)
    )

    print("\nUnique rejected lipid IDs by dataset:")
    print(rejected_counts.to_string())

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(
        rejected_counts.index,
        rejected_counts.values,
        edgecolor="black",
        color=bar_color
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Unique Rejected Lipid IDs")
    ax.set_title(
        "Lipid IDs Rejected by Fragment Validation"
    )

    for index, value in enumerate(
        rejected_counts.values
    ):
        ax.text(
            index,
            value,
            str(int(value)),
            ha="center",
            va="bottom"
        )

    ax.tick_params(
        axis="x",
        rotation=45
    )

    fig.tight_layout()

    png_file = f"{output_basename}.png"
    pdf_file = f"{output_basename}.pdf"

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    print("\nRejected-ID plot saved as:")
    print(os.path.abspath(png_file))
    print(os.path.abspath(pdf_file))

    return rejected_counts

