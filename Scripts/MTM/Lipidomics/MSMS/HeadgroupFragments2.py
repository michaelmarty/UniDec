import os
import time

import pandas as pd
import numpy as np


def str_to_float_list(s):
    try:
        return np.array([float(x) for x in s.split(" ") if x.strip() != ""])
    except:
        return np.array([])

def get_fragments_df(csvfile, allowed_classes=None):
    starttime = time.perf_counter()
    df = pd.read_csv(csvfile, dtype={"T1": str, "T2": str, "T3": str, "T4": str})

    if allowed_classes is not None:
        df = df[df["CompoundClass"].isin(allowed_classes)].reset_index(drop=True)

    # Remove [ and ] if present
    df["MSMS_mz"] = df["MSMS_mz"].str.replace("[", "").str.replace("]", "")
    df["NL_mz"] = df["NL_mz"].str.replace("[", "").str.replace("]", "")
    # Drop front and back whitespace
    df["MSMS_mz"] = df["MSMS_mz"].str.strip()
    df["NL_mz"] = df["NL_mz"].str.strip()
    # Split MSMS_mz and NL_mz columns into lists of floats
    # Use try-except to handle empty strings
    df["MSMS_mz"] = df["MSMS_mz"].apply(str_to_float_list)
    df["NL_mz"] = df["NL_mz"].apply(str_to_float_list)


    print("Created fragments dataframe:", time.perf_counter()-starttime)
    # Find all common headgroup fragments for each class and adduct
    tuples = list(zip(df["CompoundClass"], df["Adduct"], df["IonMode"]))

    unique_tuples = set(tuples)
    # Sort by class, adduct, ionmode
    unique_tuples = sorted(unique_tuples, key=lambda x: (x[0], x[1], x[2]))
    headgroup_fragments = {key: np.array([]) for key in unique_tuples}
    headgroup_fragments_nl = headgroup_fragments
    print(len(headgroup_fragments.keys()))
    print(time.perf_counter()-starttime)

    for key in headgroup_fragments.keys():
        class_name, adduct, ionmode = key
        subset = df[(df["CompoundClass"] == class_name) & (df["Adduct"] == adduct) & (df["IonMode"] == ionmode)]
        if len(subset) < 2:
            continue
        slen = len(subset)
        cutoff = slen * 0.9


        msmslist = np.concatenate(subset["MSMS_mz"].to_list())
        fragment_ions = np.round(np.ravel(msmslist), 2)
        unique, counts = np.unique(fragment_ions, return_counts=True)
        common_ions = unique[counts >= cutoff]
        headgroup_fragments[key] = common_ions

        # Repeat for neutral losses
        nllist = np.concatenate(subset["NL_mz"].to_list())
        nl_ions = np.round(np.ravel(nllist), 2)
        unique_nl, counts_nl = np.unique(nl_ions, return_counts=True)
        common_nl_ions = unique_nl[counts_nl >= cutoff]
        headgroup_fragments_nl[key] = common_nl_ions

        print(f"Class: {class_name}, Adduct: {adduct}, IonMode: {ionmode}, Common Fragments: {common_ions}, Common NLs: {common_nl_ions}")

    # print(headgroup_fragments)
    # Write to CSV
    hgfrag_list = []
    for key, fragments in headgroup_fragments.items():
        nl = headgroup_fragments_nl[key]
        if len(fragments) == 0 and len(nl) == 0:
            continue
        class_name, adduct, ionmode = key
        fstring = ",".join([str(f) for f in fragments])
        nlstring = ",".join([str(f) for f in headgroup_fragments_nl[key]])
        hgfrag_list.append({
            "CompoundClass": class_name,
            "Adduct": adduct,
            "IonMode": ionmode,
            "HeadgroupFragment_mz": fstring,
            "HeadgroupNL_mz": nlstring
        })

    hgfrag_df = pd.DataFrame(hgfrag_list)
    print(hgfrag_df)
    return hgfrag_df, df

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

poslib = "PosFullMSMS.csv"
neglib = "NegFullMSMS.csv"

negfrag, negfulldf = get_fragments_df(neglib)
posfrag, posfulldf = get_fragments_df(poslib)


# Merge posfrag and negfrag dataframes
hgfrag_df = pd.concat([posfrag, negfrag], ignore_index=True)
# Sorty by CompoundClass, Adduct, IonMode
hgfrag_df = hgfrag_df.sort_values(by=["CompoundClass", "Adduct", "IonMode"]).reset_index(drop=True)
print(hgfrag_df)

# hgfrag_df.to_csv("LipidHeadgroupFragments.csv", index=False)
# posfulldf.to_csv("PosFullMSMS.csv", index=False)
# negfulldf.to_csv("NegFullMSMS.csv", index=False)
# Convert sets to sorted lists


