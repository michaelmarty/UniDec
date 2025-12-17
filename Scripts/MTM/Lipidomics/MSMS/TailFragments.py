import os
import time

import pandas as pd
import numpy as np
from Scripts.MTM.Lipidomics.LipidFunctions import set_basic_tail_names


def str_to_float_list(s, mzlimit=0):
    try:
        a = np.array([float(x) for x in s.split(" ") if x.strip() != ""])
        a = a[a >= mzlimit]
        return a
    except:
        return np.array([])

def str_to_float_list_limit(s, mzlimit=50):
    try:
        a = np.array([float(x) for x in s.split(" ") if x.strip() != ""])
        a = a[a >= mzlimit]
        return a
    except:
        return np.array([])


def get_fragments_df(dffile):
    starttime = time.perf_counter()
    df = pd.read_csv(dffile, dtype={"T1": str, "T2": str, "T3": str, "T4": str})
    df = set_basic_tail_names(df, columnname="Name")
    df.to_csv(dffile, index=False)
    print(df.columns)

    # Remove [ and ] if present
    df["MSMS_mz"] = df["MSMS_mz"].str.replace("[", "").str.replace("]", "")
    df["NL_mz"] = df["NL_mz"].str.replace("[", "").str.replace("]", "")
    # Drop front and back whitespace
    df["MSMS_mz"] = df["MSMS_mz"].str.strip()
    df["NL_mz"] = df["NL_mz"].str.strip()
    # Split MSMS_mz and NL_mz columns into lists of floats
    # Use try-except to handle empty strings
    df["MSMS_mz"] = df["MSMS_mz"].apply(str_to_float_list_limit)
    df["NL_mz"] = df["NL_mz"].apply(str_to_float_list)

    # return df
    print("Loaded dataframe:", dffile, time.perf_counter() - starttime)
    # Find all common headgroup fragments for each class and adduct
    headgroup_fragments = []
    for index, row in df.iterrows():
        key = (row["CompoundClass"], row["Adduct"], row["IonMode"])
        if key not in headgroup_fragments:
            headgroup_fragments.append(key)
    print(len(headgroup_fragments))
    print(time.perf_counter() - starttime)

    # Find all common tails in T1, T2, T3, T4 columns, which are strings
    tails = df[["T1", "T2", "T3", "T4"]].values.ravel()
    tails = tails.astype(str)
    tails = tails[tails != "nan"]
    common_tails = np.unique(tails)
    print("Common tails:", common_tails)
    print(len(common_tails))

    prod_array = np.empty((len(headgroup_fragments), len(common_tails)), dtype=object)
    prod_array[:] = ""
    prod_array_nl = np.empty((len(headgroup_fragments), len(common_tails)), dtype=object)
    prod_array_nl[:] = ""

    hg_array = np.empty((len(headgroup_fragments),), dtype=object)
    hg_array[:] = ""
    hg_array_nl = np.empty((len(headgroup_fragments),), dtype=object)
    hg_array_nl[:] = ""

    for i, key in enumerate(headgroup_fragments):
        class_name, adduct, ionmode = key
        subset = df[(df["CompoundClass"] == class_name) & (df["Adduct"] == adduct) & (df["IonMode"] == ionmode)]
        if len(subset) < 2:
            continue

        fragment_sets = [set(np.round(mz_list, 2)) for mz_list in subset["MSMS_mz"]]
        hg_fragments = set.intersection(*fragment_sets) if fragment_sets else set()
        hg_list = np.array(sorted(hg_fragments))
        hg_string = ",".join([str(f) for f in hg_list])
        hg_array[i] = hg_string

        # Repeat for neutral losses
        nl_fragment_sets = [set(np.round(nl_list, 2)) for nl_list in subset["NL_mz"]]
        hg_nl_fragments = set.intersection(*nl_fragment_sets) if nl_fragment_sets else set()
        nl_list = np.array(sorted(hg_nl_fragments))
        hg_nl_string = ",".join([str(f) for f in nl_list])
        hg_array_nl[i] = hg_nl_string

        for j, tail in enumerate(common_tails):
            tail_subset = subset[
                (subset["T1"] == tail) |
                (subset["T2"] == tail) |
                (subset["T3"] == tail) |
                (subset["T4"] == tail)
                ]
            if len(tail_subset) < 5:
                continue

            if np.allclose(subset["MZ"], subset['MZ'].values[0], rtol=0.01):
                continue

            fragment_sets = [set(np.round(mz_list, 2)) for mz_list in tail_subset["MSMS_mz"]]
            common_fragments = set.intersection(*fragment_sets) if fragment_sets else set()
            # Make sure they are not already in headgroup fragments
            common_fragments = common_fragments - hg_fragments
            common_fragments = np.array(sorted(common_fragments))
            fstring = ",".join([str(f) for f in common_fragments])
            prod_array[i, j] = fstring

            # Repeat for neutral losses
            nl_fragment_sets = [set(np.round(nl_list, 2)) for nl_list in tail_subset["NL_mz"]]
            common_nl_fragments = set.intersection(*nl_fragment_sets) if nl_fragment_sets else set()
            # Make sure they are not already in headgroup nl fragments
            common_nl_fragments = common_nl_fragments - hg_nl_fragments
            common_nl_fragments = np.array(sorted(common_nl_fragments))
            nlstring = ",".join([str(f) for f in common_nl_fragments])
            prod_array_nl[i, j] = nlstring

            print(
                f"Class: {class_name}, Adduct: {adduct}, IonMode: {ionmode}, Tail: {tail}, Common Fragments: {common_fragments}, Common NLs: {common_nl_fragments}, HG Fragments: {hg_string}, HG NLs: {hg_nl_string}")

    # Write to CSV
    hgfrag_list = []
    for i, headgroup in enumerate(headgroup_fragments):
        class_name, adduct, ionmode = headgroup
        for j, tail in enumerate(common_tails):
            fragments = prod_array[i, j]
            nl = prod_array_nl[i, j]
            hg = hg_array[i]
            hg_nl = hg_array_nl[i]
            if fragments == "" and nl == "" and hg == "" and hg_nl == "":
                continue
            hgfrag_list.append({
                "CompoundClass": class_name,
                "Adduct": adduct,
                "IonMode": ionmode,
                "Tail": tail,
                "TailFragment_mz": fragments,
                "TailNL_mz": nl,
                "HeadgroupFragment_mz": hg,
                "HeadgroupNL_mz": hg_nl
            })

    hgfrag_df = pd.DataFrame(hgfrag_list)
    print(hgfrag_df)
    return hgfrag_df, df


os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

poslib = "PosFullMSMS.csv"
neglib = "NegFullMSMS.csv"

negfrag, negdf = get_fragments_df(neglib)
posfrag, posdf = get_fragments_df(poslib)

#
# Merge posfrag and negfrag dataframes
hgfrag_df = pd.concat([posfrag, negfrag], ignore_index=True)
# Sorty by CompoundClass, Adduct, IonMode
hgfrag_df = hgfrag_df.sort_values(by=["CompoundClass", "Adduct", "IonMode", "Tail"]).reset_index(drop=True)
print(hgfrag_df)

hgfrag_df.to_csv("LipidTailFragments.csv", index=False)
# # Convert sets to sorted lists
