from TailFragmentsAssigner import *
import time
from Scripts.MTM.Lipidomics.LipidFunctions import check_match

def assign_df(df, hgdf):
    starttime = time.perf_counter()
    outrows = []
    for index, row in df.iterrows():
        row = assign_msms(row, tol=0.05, hgdf=hgdf)
        outrows.append(row)
        if (index + 1) % 10000 == 0:
            print("Processed", index + 1, "rows of ", len(df))
            print("\tTime elapsed:", time.perf_counter() - starttime, " with estimated total time:", (time.perf_counter() - starttime) / (index + 1) * len(df))
    return pd.DataFrame(outrows)


def assign_msms(row, tol=0.05, hgdf=None, verbose=False):
    # Unpack row
    name = row["Name"]
    cname = row["CompoundClass"]
    adduct = row["Adduct"]
    tail1 = row["T1"]
    tail2 = row["T2"]
    tail3 = row["T3"]
    tail4 = row["T4"]
    tails = [tail1, tail2, tail3, tail4]
    msms_mz = row["MSMS_mz"]
    nl_mz = row["NL_mz"]

    # Get these from the assignment module
    ignorefrag = ignore_fragments[cname] if cname in ignore_fragments else []
    ignorenl = ignore_nls[cname] if cname in ignore_nls else []

    fragdict, nldict = gen_fragdict(cname, adduct, tails, hgdf)

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
