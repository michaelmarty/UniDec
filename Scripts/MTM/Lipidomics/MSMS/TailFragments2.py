import os
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
print(mass_slfrag)


# CER ref: https://link.springer.com/article/10.1007/s00216-013-7601-y

hgbasedict = {"PG": [228.04], "PI": [316.05], "SQDG": [282.04], "PS": [87.032, 223.0022], "OxPS": [87.034, 223.0022],
              "PC": [281.102], "OxPG": [228.04], "MGDG": [218.077], "EtherPG": [228.04 - mass_o],
              "EtherPI": [316.05 - mass_o, 316.05, 162.05], "EtherPS": [87.032 - mass_o, 223.0022 - mass_o, 87.032],
              "DGGA": [250.07], "Ac2PIM1": [460.1, 136.00], "DLCL": [210.03, 364.032], "AHexCer": [208.052, 222.08]}

ignore_fragments = {}
ignore_nls = {"Hex2Cer": [46.01], "Hex3Cer": [46.01], "HexCer": [46.01], "AHexCer": [46.01]}
acelated_classes = ["AHexCer"]



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
        if numox == 2:
            mass += mass_n + 3 * mass_h
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
        fragments["FA(+O-H2O+CH2)"] = mass - mass_h2o + mass_ch2
        fragments["FA(-H)"] = acylmass - mass_h
        fragments["FA(-H2O+CH2)"] = acylmass -mass_h - mass_h2o + mass_ch2
        fragments["FA(-CH2)"] = acylmass - mass_ch2 - mass_h

        if nhydrox > 0:
            fragments["FA(-CH2O)"] = acylmass - mass_ch2 - mass_o
            fragments["FA(+O-2H)"] = mass - 2 * mass_h


        # fragments["FA+Glycerol"] = mass + mass_glycerol - mass_h2o  # FA + glycerol - H2O
        fragments[
            "FA+GP"] = mass + mass_glycerol + mass_phosphate - 2 * mass_h2o + 2 * mass_h  # FA + glycerol + phosphate - 2 H2O + 2 H
        fragments[
            "FA+GP(-H2O)"] = mass + mass_glycerol + mass_phosphate - 3 * mass_h2o + 2 * mass_h  # FA + glycerol + phosphate - 3 H2O

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

        if nhydrox > 0:
            nls["FA(+O+C2H3O2)"] = mass + 2 * mass_c + 2 * mass_o + 3 * mass_h
            nls["FA(+O+C2H7O2)"] = mass + 2 * mass_c + 2 * mass_o + 7 * mass_h
            nls["FA(+O+C2H5O3)"] = mass + 2 * mass_c + 3 * mass_o + 5 * mass_h
            nls["FA(+O+C2H9O3)"] = mass + 2 * mass_c + 3 * mass_o + 9 * mass_h
            nls["FA(+O+C3H7O3)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h
            nls["FA(+O+C4H10O3N1)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h
            nls["FA(+O+C2H3O2-CH2)"] = mass + 2 * mass_c + 2 * mass_o + 3 * mass_h - mass_ch2
            nls["FA(+O+C2H7O2-CH2)"] = mass + 2 * mass_c + 2 * mass_o + 7 * mass_h- mass_ch2
            nls["FA(+O+C2H5O3-CH2)"] = mass + 2 * mass_c + 3 * mass_o + 5 * mass_h- mass_ch2
            nls["FA(+O+C2H9O3-CH2)"] = mass + 2 * mass_c + 3 * mass_o + 9 * mass_h- mass_ch2
            nls["FA(+O+C3H7O3-CH2)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h- mass_ch2
            nls["FA(+O+C4H10O3N1-CH2)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h- mass_ch2
            nls["FA(+O-H)"] = mass - mass_h
            nls["FA(+O+3H)"] = mass +3 * mass_h
            nls["FA(+O-H+H2O)"] = mass - mass_h + mass_h2o
            nls["FA(+O+3H+H2O)"] = mass + 3 * mass_h + mass_h2o
            nls["FA(+O+CH3O)"] = mass + mass_c + 3 * mass_h + mass_o
            nls["FA(+O+C2H6ON)"] = mass + 2 * mass_c + mass_n + 6 * mass_h + mass_o

            # nls["FA(+O+C3H7O3-CH2)"] = mass + 3 * mass_c + 3 * mass_o + 7 * mass_h
            # nls["FA(+O+C4H10O3N1-CH2)"] = mass + 4 * mass_c + 3 * mass_o + mass_n + 10 * mass_h - mass_ch2


        if adduct == "[M+CH3COO]-":
            # For Ceramides
            nls["FA(+O+H2O)+Acetate"] = mass + mass_h2o + mass_acetate + 1 * mass_h
            nls[
                "FA(+NH+C2H2+H3O)+Acetate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_h2o + mass_h + + mass_acetate + 1 * mass_h
            nls[
                "FA(+NH+C3H4O+H3O)+Acetate"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o + mass_h2o + mass_h + mass_acetate + 1 * mass_h

        if adduct == "[M+HCOO]-":
            nls["FA(+O+H2O)+Formate"] = mass + mass_h2o + mass_formate + 1 * mass_h
            nls[
                "FA(+NH+C2H2+H3O)+Formate"] = mass - mass_o + mass_n + 3 * mass_h + 2 * mass_c + mass_h2o + mass_h + mass_formate + 1 * mass_h
            nls[
                "FA(+NH+C3H4O+H3O)+Formate"] = mass - mass_o + mass_n + 5 * mass_h + 3 * mass_c + mass_o + mass_h2o + mass_h + mass_formate + 1 * mass_h

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

        outstring2 = ",".join(outnllist)

    row["TailFragment_Assigned"] = outstring
    row["TailNL_Assigned"] = outstring2

    # if "NOTFOUND" in outstring or "NOTFOUND" in outstring2:
    #     print(row)

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


# print(predict_tail_fragments("18:0(methyl)"))
# exit()

if __name__ == "__main__":

    import time

    starttime = time.perf_counter()
    os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

    df = pd.read_csv("LipidTailFragments.csv", dtype=str)

    # exit()

    # Select df for negative ion mode
    df = df[df["IonMode"] == "negative"].reset_index(drop=True)
    # Select df for specific class
    # df = df[df["CompoundClass"] == "Cer_ADS"].reset_index(drop=True)
    # Isolate tails with 18:1
    # df_18_1 = df[df["Tail"] == "16:0"].reset_index(drop=True)

    df = assign_fragments_to_df(df, tol=0.05)
    df = merge_fragments_df(df)
    # Drop rows where both TailFragment_Assigned and TailNL_Assigned are empty
    # df = df[
    #     ~((df["TailFragment_Assigned"].isna() | (df["TailFragment_Assigned"].str.strip() == "")))].reset_index(
    #     drop=True)
    print(df)
    # Write to xlsx
    df.to_excel("LipidTailFragments_Negative.xlsx", index=False)
    print("Written output file.")
    # Send df to list and print

    print(time.perf_counter() - starttime)
