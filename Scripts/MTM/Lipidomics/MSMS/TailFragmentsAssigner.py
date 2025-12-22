import os
import pandas as pd
import styleframe as sf
from PublicScripts.Lipids.LipidFragPrediction import ignore_nls, predict_tail_fragments

def assign_fragments(row, tol=0.05):
    tail = row["Tail"]
    fraglist = row["TailFragment_mz"]
    nllist = row["TailNL_mz"]
    hgfrags = row["HeadgroupFragment_mz"]
    cname = row["CompoundClass"]
    adduct = row["Adduct"]
    mode = row["IonMode"]
    ignorenl = ignore_nls[cname] if cname in ignore_nls else []
    # print(cname, adduct)

    predicted_frags, predicted_nls = predict_tail_fragments(tail, hgfrags, adduct=adduct, classname=cname,mode=mode)
    outfraglist = []
    outnllist = []

    if not (pd.isna(fraglist) or fraglist.strip() == ""):
        fraglist = [float(f) for f in fraglist.split(",")]
        for frag in fraglist:
            found = False
            for frag_name, frag_mz in predicted_frags.items():
                if found:
                    continue
                if abs(frag - frag_mz) <= tol:
                    outfraglist.append(frag_name)
                    found = True
            if not found:
                outfraglist.append("NOTFOUND(" + str(frag) + "|" + tail + ")")
                # print("NOTFOUND Fragment:", frag, "for", cname, "Tail:", tail)
                # defect = frag - predicted_frags['FA(+OH)']
                # print(defect)
                # form = calc_fromula_from_mass([defect], searchspace=reallysmallsearch, tolerance=200)
                # print(form)
                # print(frag - predicted_frags['FA(+OH)'], ",")

    if not (pd.isna(nllist) or nllist.strip() == ""):
        nllist = [float(nl) for nl in nllist.split(",")]
        for nl in nllist:
            if nl in ignorenl:
                outnllist.append("IGNORE(" + str(nl) + "|" + tail + ")")
                continue
            found = False
            for nl_name, nl_mz in predicted_nls.items():
                if found:
                    continue
                if abs(nl - nl_mz) <= tol:
                    outnllist.append(nl_name)
                    found = True
            if not found:
                outnllist.append("NOTFOUND(" + str(nl) + "|" + tail + ")")
                # print("NOTFOUND Neutral Loss:", nl, "for", cname, "Tail:", tail)
                # defect = nl - predicted_nls['FA(+OH)']
                # print(defect)
                # form = calc_fromula_from_mass([defect], searchspace=reallysmallsearch, tolerance=200)
                # print(form)
                # print(nl - predicted_frags['FA(+OH)'] + mass_h, ",")

    if len(outfraglist) == 0 and len(outnllist) == 0:
        outstring = ""
        outstring2 = ""
    else:
        # outfraglist, outnllist = merge_lists(outfraglist, outnllist)
        outstring = ",".join(outfraglist)
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
    df = df[df["IonMode"] == "positive"].reset_index(drop=True)
    # Select df for specific class
    # df = df[df["CompoundClass"] == "SL"].reset_index(drop=True)
    # Isolate tails with 18:1
    # df = df[df["Tail"] == "18:1;2O"].reset_index(drop=True)

    df = assign_fragments_to_df(df, tol=0.05)
    df = merge_fragments_df(df)
    # Drop rows where both TailFragment_Assigned and TailNL_Assigned are empty
    # df = df[
    #     ~((df["TailFragment_Assigned"].isna() | (df["TailFragment_Assigned"].str.strip() == "")))].reset_index(
    #     drop=True)
    # print(df)
    # Write to xlsx
    # df.to_excel("LipidTailFragments_18120.xlsx", index=False, na_rep='NaN')
    # print("Written output file.")
    # Send df to list and print

    writer = sf.ExcelWriter("LipidTailFragments_Positive.xlsx")
    sfdf = sf.StyleFrame(df)
    # df.to_excel(writer, sheet_name='Sheet 1', index=False, na_rep='NaN')
    sfdf.to_excel(writer, best_fit=df.columns.tolist())

    writer.close()

    print(time.perf_counter() - starttime)
