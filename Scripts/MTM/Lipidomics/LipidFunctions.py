import numpy as np
import re

adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M+HCOO]-", "[M+HCOO]", "-1", -1, "[M+HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"],
                              ["[M+CH3COO]-", "[M+CH3COO]", "-1", -1, "[M+CH3COO]1-"],
                              ["[M+Na]+", "[M+Na]", "1", -1, "[M+Na]1+"]])

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


def lipidlist_to_LC(df):
    outnames = []
    fixednames = []
    fixedclasses = []
    for i, row in df.iterrows():
        lclass = row["Molecule List Name"]
        name = row["Molecule Name"]
        name = fix_dash(name)

        name, lclass = fix_name(name, lclass)

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
    if "TG" in name:
        name = name.replace("TG", "TAG")
        lclass = "TAG"

    if "DG" in name:
        name = name.replace("DG", "DAG")
        lclass = "DAG"

    if "MG" in name:
        name = name.replace("MG", "MAG")
        lclass = "MAG"

    return name, lclass


def fix_d7(name, adduct):
    if "d7" in name:
        adduct = adduct.replace("M", "M7H2")
        name= name.replace("(d7)", " + D7")

    if "d9" in name:
        adduct = adduct.replace("M", "M9H2")
        name = name.replace("(d9)", " + D9")

    return name, adduct


def count_tails(name, verbose=False):
    # Remove values of just a name or a name and then an adduct
    if " " not in name or ":" not in name:
        return "NA", -1, -1

    try:
        tname = name.split(" ")[1]
    except Exception as e:
        print("Unable to parse tails: ", name)
        raise e

    # fas = tname.split("_")
    fas = re.split('_|/', tname)
    tls = 0
    tus = 0
    # Remove '' from fas
    for f in fas:
        if f == '':
            continue
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
    return tname, tl, tu


def set_tails(df):
    tails = []
    lengths = []
    unsats = []
    for i, row in df.iterrows():
        name = row["Molecule Name"]
        tname, tl, tu = count_tails(name)
        tails.append(tname)
        lengths.append(tl)
        unsats.append(tu)
    df["Tails"] = tails
    df["Tail Length"] = lengths
    df["Tail Unsaturation"] = unsats
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
