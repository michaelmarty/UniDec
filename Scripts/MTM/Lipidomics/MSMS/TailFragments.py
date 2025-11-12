import os
import time

import pandas as pd
import numpy as np
from Scripts.MTM.Lipidomics.LipidFunctions import set_basic_tail_names

def get_fragments_df(dffile):
    starttime = time.perf_counter()
    df = pd.read_csv(dffile, dtype={"T1": str, "T2": str, "T3": str, "T4": str})
    # df = set_basic_tail_names(df, columnname="Name")
    # df.to_csv(dffile, index=False)
    # return df
    print("Loaded dataframe:", dffile, time.perf_counter()-starttime)
    # Find all common headgroup fragments for each class and adduct
    headgroup_fragments = {}
    headgroup_fragments_nl = {}
    for index, row in df.iterrows():
        key = (row["CompoundClass"], row["Adduct"], row["IonMode"])
        if key not in headgroup_fragments:
            headgroup_fragments[key] = set()
            headgroup_fragments_nl[key] = set()
    print(headgroup_fragments.keys())
    print(time.perf_counter()-starttime)

    # Find all common tails in T1, T2, T3, T4 columns, which are strings
    common_tails = np.unique(df[["T1", "T2", "T3", "T4"]].values.ravel())
    print("Common tails:", common_tails)

    exit()

    for key in headgroup_fragments.keys():
        class_name, adduct, ionmode = key
        subset = df[(df["CompoundClass"] == class_name) & (df["Adduct"] == adduct) & (df["IonMode"] == ionmode)]
        if len(subset) < 2:
            continue
        fragment_sets = [set(np.round(mz_list, 2)) for mz_list in subset["MSMS_mz"]]
        common_fragments = set.intersection(*fragment_sets) if fragment_sets else set()
        common_fragments = np.array(sorted(common_fragments))
        headgroup_fragments[key] = common_fragments

        # Repeat for neutral losses
        nl_fragment_sets = [set(np.round(nl_list, 2)) for nl_list in subset["NL_mz"]]
        common_nl_fragments = set.intersection(*nl_fragment_sets) if nl_fragment_sets else set()
        common_nl_fragments = np.array(sorted(common_nl_fragments))
        headgroup_fragments_nl[key] = common_nl_fragments

        print(f"Class: {class_name}, Adduct: {adduct}, IonMode: {ionmode}, Common Fragments: {common_fragments}, Common NLs: {common_nl_fragments}")

    print(headgroup_fragments)
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

posfrag = get_fragments_df(poslib)
negfrag = get_fragments_df(neglib)
#
# # Merge posfrag and negfrag dataframes
# hgfrag_df = pd.concat([posfrag, negfrag], ignore_index=True)
# # Sorty by CompoundClass, Adduct, IonMode
# hgfrag_df = hgfrag_df.sort_values(by=["CompoundClass", "Adduct", "IonMode"]).reset_index(drop=True)
# print(hgfrag_df)
#
# hgfrag_df.to_csv("LipidTailFragments.csv", index=False)
# # Convert sets to sorted lists


