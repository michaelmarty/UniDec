import os
import pandas as pd
from Scripts.MTM.Lipidomics.LipidFunctions import set_basic_tail_names, check_match
from TailFragmentsAssigner import *

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
    matchednl = ref_mz - matchedmz

    # Filter both to >5 m/z
    matchednl = matchednl[matchednl > 5]
    matchedmz = matchedmz[matchedmz < ref_mz-5]

    # Get the potential adduct neutral losses
    adductnls = adduct_nls[cname] if cname in adduct_nls else []

    fragdict, nldict = gen_fragdict(cname, adduct, tails, hgdf)

    if verbose:
        print(name, fragdict, nldict)


    nfrag_missed = 0
    nnl_missed = 0
    tailcounts = np.zeros(4, dtype=int)
    hgcounts = 0
    adductcounts = 0

    outfraglist = []
    for frag in matchedmz:
        found = False
        for frag_name, frag_mz in fragdict.items():
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
            if verbose:
                print("NOTFOUND Fragment:", frag, "for", cname, "Name:", name)
    outstring = ",".join(outfraglist)

    tailnlcounts = np.zeros(4, dtype=int)
    hgcountsnl = 0
    adductcountsnl = 0
    outnllist = []
    for nl in matchednl:
        if nl in adductnls:
            outnllist.append("Adduct[" + str(nl) + "]")
            adductcountsnl += 1
            continue
        found = False
        for nl_name, nl_mz in nldict.items():
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
            if verbose:
                print("NOTFOUND Neutral Loss:", nl, "for", cname, "Name:", name)
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

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

negfile = "PosIDsC1_Ref.csv"

negdf = pd.read_csv(negfile)

hgfile = "LipidHeadgroupFragments.csv"
hgdf = pd.read_csv(hgfile, dtype=str)

outdf = assign_df(negdf, hgdf=hgdf)

outdf.to_csv(negfile[:-4] + "_TailCheck.csv", index=False)