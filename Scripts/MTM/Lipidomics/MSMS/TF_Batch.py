from TailFragmentsAssigner import *
import time

def assign_df(df, hgdf):
    starttime = time.perf_counter()
    outrows = []
    for index, row in df.iterrows():
        row = assign_msms(row, tol=0.05, hgdf=hgdf)
        outrows.append(row)
        if (index + 1) % 10000 == 0:
            print("Processed", index + 1, "rows of ", len(df))
            print("\tTime elapsed:", time.perf_counter() - starttime, " with estimated total time:", (time.perf_counter() - starttime) / (index + 1) * len(df))

def check_match(query_mz, target, tol=0.05):
    if abs(query_mz - target) <= tol:
        return True
    return False

def assign_msms(row, tol=0.05, hgdf=None, verbose=False):
    # Unpack row
    name = row["Name"]
    cname = row["CompoundClass"]
    adduct = row["Adduct"]
    tail1 = row["T1"]
    tail2 = row["T2"]
    tail3 = row["T3"]
    tail4 = row["T4"]
    msms_mz = row["MSMS_mz"]
    nl_mz = row["NL_mz"]

    # Get these from the assignment module
    ignorefrag = ignore_fragments[cname] if cname in ignore_fragments else []
    ignorenl = ignore_nls[cname] if cname in ignore_nls else []
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
    tails = [tail1, tail2, tail3, tail4]
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

    if verbose:
        print(name, fragdict, nldict)

    nfrag_missed = 0
    nnl_missed = 0
    if pd.isna(msms_mz) or msms_mz.strip() == "":
        outstring = ""
    else:
        msms_mz = msms_mz.replace("[", "").replace("]", "")
        mzsplit = msms_mz.split(" ")
        mzsplit = [f for f in mzsplit if f.strip() != ""]
        fraglist = [float(f) for f in mzsplit]

        outfraglist = []
        for frag in fraglist:
            if frag in ignorefrag:
                continue

            found = False
            for frag_name, frag_mz in fragdict.items():
                if check_match(frag, frag_mz, tol):
                    outfraglist.append(frag_name)
                    found = True
            if not found:
                outfraglist.append("NOTFOUND[" + str(frag) + "|" + name + "]")
                nfrag_missed += 1
                if verbose:
                    print("NOTFOUND Fragment:", frag, "for", cname, "Name:", name)

        outstring = ",".join(outfraglist)

    if pd.isna(nl_mz) or nl_mz.strip() == "":
        outstring2 = ""
    else:
        nl_mz = nl_mz.replace("[", "").replace("]", "")
        nlsplit = nl_mz.split(" ")
        nlsplit = [nl for nl in nlsplit if nl.strip() != ""]
        nllist = [float(nl) for nl in nlsplit]

        outnllist = []
        for nl in nllist:
            if nl in ignorenl:
                continue
            found = False
            for nl_name, nl_mz in nldict.items():
                if check_match(nl, nl_mz, tol):
                    outnllist.append(nl_name)
                    found = True
            if not found:
                outnllist.append("NOTFOUND[" + str(nl) + "|" + name + "]")
                nnl_missed += 1
                if verbose:
                    print("NOTFOUND Neutral Loss:", nl, "for", cname, "Name:", name)
        outstring2 = ",".join(outnllist)

    row["Fragments_Assigned"] = outstring
    row["NL_Assigned"] = outstring2
    row["Fragments_Missed"] = nfrag_missed
    row["NL_Missed"] = nnl_missed
    return row


os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

negfile = "NegFullMSMS.csv"

hgfile = "LipidHeadgroupFragments.csv"
hgdf = pd.read_csv(hgfile, dtype=str)

negdf = pd.read_csv(negfile, dtype=str)

outdf = assign_df(negdf, hgdf)

outdf.to_csv("NegFullMSMS_Assigned.csv", index=False)
