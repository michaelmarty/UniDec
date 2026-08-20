import pandas as pd
import os
from io import StringIO
from PublicScripts.Lipids.LipidFragPrediction import *
from PublicScripts.Lipids.MatchPlotting import *
from PublicScripts.Lipids.LipidFunctions import parse_spec_string, fix_adduct, match_frag_to_spec, set_tails, \
    spec_to_nl, write_colored_excel, flag_weird_lipids, get_tail_carbons, get_tail_unsaturation
import numpy as np
import re
import ms_entropy as me

# Define CSV

lipid_data = """
Class	Adduct	Required Frags	Disqualifying Frags

Cer	[M+H]	FA1(+O-H2O), FA1(+O-2H2O), FA1(+O-CH4O2)	
Cer	[M+H]	FA1(+O-H2O), FA1(+O-2H2O), FA1(+O-3H2O)	
Cer	[M-H]	FA1(+O-SL)	
Cer	[M-H]	FA1(+O-SL+C2H2)	
Cer	[M-H]	FA1(+O-SL), FA2(+O-H2O)	
Cer	[M-H]	FA1(+O-SL+C2H2), FA2(+O-H2O)	
Cer	[M+CH3COO]	FA1(+O-SL), FA2(+O-H2O)	
Cer	[M+CH3COO]	FA1(+O-SL), FA2(+O)	
CL	[M-H]	FA1(+O), FA2(+O), FA3(+O), FA4(+O), HG(153)
CL	[M+NH4]	HGNL(17.03)
CoQ	[M+H]	HG(197)	
DG	[M+NH4]	FA1(+G-OH), FA2(+G-OH)
FAHFA	[M-H]	FA1(+O)	HG(153), HG(171), HG(196), HG(241)
GM3	[M-H]	HG(290), HGNL(291)	
HexCer	[M+H]	FA1(+O-H2O), FA1(+O-CH4O2)	
HexCer	[M+H]	FA1(+O), FA1(+O-H2O), FA1(+O-3H2O)	
HexCer	[M-H]	FA1(+O-SL), FA2(+O-H2O)	
HexCer	[M-H]	FA2(+O-CH2O), FA2(+O)	
HexCer	[M+CH3COO]	FA1(+O-SL), FA2(+O-H2O)	
Hex2Cer	[M+H]	HGNL(342), FA1(+O-H2O), FA1(+O-CH4O2)	
Hex2Cer	[M+CH3COO]	HG(179), HGNL(222), HGNL(384)	
LPA	[M-H]	HG(153)	
LPC	[M+H]	HG(184)	
LPC	[M+CH3COO]	HG(224), FA1(+O)	
EtherLPC	[M+CH3COO]	HG(168)	
EtherLPC	[M+H]	HG(184)	
LPE	[M-H]	HG(196), FA1(+O)	
LPG	[M-H]	HG(153), FA1(+O)	
LPI	[M-H]	HG(153), FA1(+O), HG(241)	
LPS	[M-H]	HG(153), FA1(+O)	
PA	[M-H]	HG(153), FA1(+O), FA2(+O)	
PC	[M+H]	HG(184)	
PC	[M+CH3COO]	HG(224), HGNL(74), FA1(+O), FA2(+O)
PC	[M+Na]	NLFA1(+O-C3H9N), NLFA2(+O-C3H9N)
PC	[M+Na]	NLFA1(+O-C5H14NO4P), NLFA2(+O-C5H14NO4P)
EtherPC	[M+CH3COO]	HG(224), FA1(+O), FA2(+O), FA1(+O-C2H4), FA2(+O-C2H4)	
EtherPC	[M+H]	HG(184)	
PE	[M+H]	HGNL(141)	
PE	[M-H]	HG(196), FA1(+O), FA2(+O)	
EtherPE	[M+H]	HGNL(141)	
EtherPE	[M-H]	FA2(+O), NLFA2(-H)	
PG	[M-H]	HG(153), FA1(+O), FA2(+O)	
PI	[M-H]	HG(241), FA1(+O), FA2(+O)	
PS	[M+H]	HGNL(185)	
PS	[M-H]	HGNL(87), FA1(+O), FA2(+O)	
SM	[M+CH3COO]	HGNL(74), HG(168)	FA1(+O), FA2(+O)
SM	[M+H]	HG(184)	FA1(+O), FA2(+O)
SM	[M+Na]	HGNL(59)
SM	[M+Na]	HGNL(183)
TG	[M+NH4]	NLFA1(+O+H2O), NLFA2(+O+H2O), NLFA3(+O+H2O)
TG	[M+Na]	NLFA1(+OH), NLFA2(+OH), NLFA3(+OH)
"""
# Read and clean the fragment-rule table
lipid_df = pd.read_csv(StringIO(lipid_data), sep="\t", index_col=False)

for col in lipid_df.columns:
    if lipid_df[col].dtype == "object":
        lipid_df[col] = lipid_df[col].str.strip()


def calc_spectral_entropy(spectrum):
    entropy = me.calculate_spectral_entropy(spectrum, clean_spectrum=False)
    return entropy


def spectral_entropy_score(spec1, spec2, tol=0.7):
    # unweighted_similarity = me.calculate_unweighted_entropy_similarity(np.array(spec1), np.array(spec2), clean_spectra=False)
    # print(f"Unweighted entropy similarity: {unweighted_similarity}.")

    # Calculate entropy similarity.
    similarity = me.calculate_entropy_similarity(np.array(spec1), np.array(spec2), ms2_tolerance_in_da=tol)
    return similarity


def get_hg_frag(lipid_class, adduct, fragname):
    # Extract the nominal fragment label, e.g. HG(184) -> "184"
    if "(" not in fragname or ")" not in fragname:
        print(f"No fragment number specified in {fragname}.")
        return None

    fragnum = fragname.split("(")[-1].split(")")[0].strip()

    # Find the correct class/adduct entry
    row = hgdf_default[
        (hgdf_default["CompoundClass"] == lipid_class)
        & (
                hgdf_default["Adduct"]
                == fix_adduct(adduct, 1, 0)
        )
        ]

    if not row.empty:
        if "HGNL" in fragname:
            vals = row["HeadgroupNL_mz"].values[0]
        elif "HG" in fragname:
            vals = row["HeadgroupFragment_mz"].values[0]
        else:
            raise ValueError(
                f"Fragment name {fragname} not recognized."
            )

        if pd.notna(vals) and str(vals).strip() != "":
            for value in str(vals).split(","):
                value = value.strip()

                # Match nominal label to exact-mass entry
                if value.startswith(fragnum):
                    return float(value)

    # Fallback only when no database value exists
    try:
        return float(fragnum)
    except ValueError:
        print(
            f"No fragment mass found for "
            f"{lipid_class}, {adduct}, {fragname}."
        )
        return None


def get_tail_frag(row, fragname):
    nlmode = False
    if fragname[:2] == "NL":
        nlmode = True
        fragname = fragname[2:]

    match = re.search(r"FA(\d+)", fragname)

    if match is None:
        return None

    index = int(match.group(1)) - 1
    tail_col = f"T{index + 1}"

    if tail_col not in row.keys():
        print(f"Tail column {tail_col} not found in row for {row['Metabolite name']}.")
        return None

    tail = row[tail_col]

    print("DEBUG")
    print("Metabolite:", row["Metabolite name"])
    print("Required frag:", fragname)
    print("Tail col:", tail_col)
    print("Tail value:", repr(tail))

    if pd.isna(tail):
        print(f"No tail information found for {row['Metabolite name']} at {tail_col}.")
        return None

    adduct = fix_adduct(row["Adduct type"], 1, 0)
    lclass = row["Ontology"]

    if adduct[-1] == "+":
        polarity = "positive"
    elif adduct[-1] == "-":
        polarity = "negative"
    else:
        print(f"Adduct {adduct} does not have a recognized polarity for {row['Metabolite name']}.")
        return None

    tail_frags, tail_nls = predict_tail_fragments(
        tail,
        adduct=adduct,
        classname=lclass,
        mode=polarity
    )

    fragname_key = "FA" + fragname[3:]

    if nlmode:
        tdict = tail_nls
    else:
        tdict = tail_frags

    print("fragname_key:", fragname_key)
    print("Available keys:", list(tdict.keys()))

    if fragname_key in tdict.keys():
        return tdict[fragname_key]

    print(f"Fragment key {fragname_key} not found for {row['Metabolite name']}")
    print(f"Tail used: {tail_col} = {tail}")
    return None

    adduct = fix_adduct(row["Adduct type"], 1, 0)
    lclass = row["Ontology"]

    if adduct[-1] == "+":
        polarity = "positive"
    elif adduct[-1] == "-":
        polarity = "negative"
    else:
        print(f"Adduct {adduct} does not have a recognized polarity for {row['Metabolite name']}.")
        return None

    tail_frags, tail_nls = predict_tail_fragments(
        tail,
        adduct=adduct,
        classname=lclass,
        mode=polarity
    )

    # Drop the number after the FA in the fragname to match keys in tail_frags
    fragname_key = "FA" + fragname[3:]
    if nlmode:
        tdict = tail_nls
    else:
        tdict = tail_frags

    if fragname_key in tdict.keys():
        frag = tdict[fragname_key]
        # print("Fragment found:", frag, tail)
        return frag


def get_match_result(spec, ref_spec, frag, tol):
    matched_frag = match_frag_to_spec(spec, frag, tol)
    matched_ref = match_frag_to_spec(ref_spec, frag, tol)

    if matched_frag and matched_ref:
        matched_correctly = "Yes"
        match_bool = True
    elif not matched_frag and matched_ref:
        matched_correctly = "No"
        match_bool = False
    elif matched_frag and not matched_ref:
        matched_correctly = "Yes but Ref Failed"
        match_bool = False
    else:
        matched_correctly = "No and Ref Failed"
        match_bool = False

    return matched_correctly, match_bool


def check_disqualifying_fragments(row, lipid_class, adduct, bad_frags, spec, ref_spec, mz, tol):
    """
    Returns:
        bad_matches_text: string describing disqualifying fragments found
        has_bad_fragment: True if any disqualifying fragment is present in the sample spectrum
    """
    if pd.isna(bad_frags) or str(bad_frags).strip() == "":
        return "", False

    bad_matches_text = ""
    has_bad_fragment = False

    for f in str(bad_frags).split(","):
        f = f.strip()
        if f == "":
            continue

        frag = None

        if "HG" in f:
            frag = get_hg_frag(lipid_class, adduct, f)

        elif "FA" in f:
            frag = get_tail_frag(row, f)

        if lipid_class == "TG" and adduct == "[M+Na]":
            print("\nSODIUM TG DEBUG")
            print("Lipid:", row["Metabolite name"])
            print("Precursor used:", mz)
            print("Fragment rule:", f)
            print("Predicted neutral loss:", frag)

            sample_nl = spec_to_nl(spec, mz)
            reference_nl = spec_to_nl(ref_spec, mz)

            print("Sample product peaks:")
            print(spec)

            print("Reference product peaks:")
            print(ref_spec)

            print("Sample neutral losses:")
            print(sample_nl)

            print("Reference neutral losses:")
            print(reference_nl)

            print(
                "Sample match:",
                match_frag_to_spec(sample_nl, frag, tol)
            )

            print(
                "Reference match:",
                match_frag_to_spec(reference_nl, frag, tol)
            )

        if frag is None:
            bad_matches_text += f"{f}: Not Found; "
            continue

        if "HGNL" in f or "NLFA" in f:
            matched_frag = match_frag_to_spec(spec_to_nl(spec, mz), frag, tol)
        else:
            matched_frag = match_frag_to_spec(spec, frag, tol)

        if matched_frag:
            bad_matches_text += f"{f}: Present - Disqualifies ID; "
            has_bad_fragment = True
        else:
            bad_matches_text += f"{f}: Not Present; "

    return bad_matches_text, has_bad_fragment


def check_fragments(df, ldf=None, plot_individual=False):
    if ldf is None:
        ldf = lipid_df

    # Loop through each row in the input DataFrame
    results = []
    bresults = []
    bad_results = []
    for index, row in df.iterrows():
        lipid_class = row["Ontology"].strip()
        adduct = row["Adduct type"]
        mz = row["Reference m/z"]
        # Clean adduct
        adduct = fix_adduct(adduct, 0, 1)

        # Match ldf to get required fragments
        required_frags = ldf[(ldf["Class"] == lipid_class) & (ldf["Adduct"] == adduct)]["Required Frags"].values
        disqualifying_frags = ldf[(ldf["Class"] == lipid_class) & (ldf["Adduct"] == adduct)][
            "Disqualifying Frags"].values

        if len(required_frags) == 0:
            print(f"No required fragments found for {lipid_class} with adduct {adduct}.")
            results.append("No required fragments found for this type")
            bresults.append(False)
            bad_results.append("")
            continue

        tol = row["Tolerance"]
        spec = parse_spec_string(row["MS/MS spectrum"], threshold=0.01, norm=True)
        ref_spec = parse_spec_string(row["Ref Spec"], threshold=0.01, norm=True)

        if len(spec) == 0 or len(ref_spec) == 0:
            print(f"Empty spectrum for {row['Metabolite name']}. Skipping.")
            results.append("Empty spectrum")
            bresults.append(False)
            bad_results.append("")
            continue

        matches = ""
        mbool = False

        for frags in required_frags:
            matches_sub = ""
            mbool_sub = True
            for f in frags.split(","):
                f = f.strip()
                if "HG" in f:
                    frag = get_hg_frag(lipid_class, adduct, f)
                    if frag is None:
                        print(f"Could not find fragment for {lipid_class} with adduct {adduct} and fragment {f}.")
                        matches_sub += f"{f}: Not Found; "
                        mbool_sub = False
                        continue
                    if "HGNL" in f:
                        matched_correctly, match_bool = get_match_result(spec_to_nl(spec, mz), spec_to_nl(ref_spec, mz),
                                                                         frag, tol)
                    else:
                        matched_correctly, match_bool = get_match_result(spec, ref_spec, frag, tol)
                    mbool_sub = mbool_sub and match_bool

                    print(f"HG fragment for {lipid_class} with adduct {adduct}: {frag}", matched_correctly)

                    matches_sub += f"{f}: {matched_correctly}; "

                if "FA" in f:
                    frag = get_tail_frag(row, f)
                    if frag is None:
                        print(f"Could not find fragment for {lipid_class} with adduct {adduct} and fragment {f}.")
                        matches_sub += f"{f}: Not Found; "
                        mbool_sub = False
                        continue

                    if "NLFA" in f:
                        matched_correctly, match_bool = get_match_result(spec_to_nl(spec, mz), spec_to_nl(ref_spec, mz),
                                                                         frag, tol)
                    else:
                        matched_correctly, match_bool = get_match_result(spec, ref_spec, frag, tol)
                    mbool_sub = mbool_sub and match_bool

                    print(f"FA fragment for {lipid_class} with adduct {adduct}: {frag}", matched_correctly)

                    matches_sub += f"{f}: {matched_correctly}; "

            if mbool_sub:
                mbool = mbool_sub
                matches = matches_sub
                continue

            if matches == "":
                matches += matches_sub

        # Check disqualifying fragments
        bad_match_text = ""
        has_bad_fragment = False

        if "Disqualifying Frags" in ldf.columns and len(disqualifying_frags) > 0:
            for bad_frags in disqualifying_frags:
                bad_sub_text, bad_sub_bool = check_disqualifying_fragments(
                    row,
                    lipid_class,
                    adduct,
                    bad_frags,
                    spec,
                    ref_spec,
                    mz,
                    tol
                )

                bad_match_text += bad_sub_text

                if bad_sub_bool:
                    has_bad_fragment = True

        if has_bad_fragment:
            mbool = False
            matches += " DISQUALIFYING FRAGMENTS FOUND: " + bad_match_text

        # print("Score based on Spectral Entropy:", spectral_entropy_score(spec, ref_spec, tol=tol))
        #
        if plot_individual and not mbool:
            butterfly_plot(spec, ref_spec, title=row["Metabolite name"])
            plt.show()

        results.append(matches)
        bresults.append(mbool)
        bad_results.append(bad_match_text)

    df["Matched Required Fragments"] = results
    df["Disqualifying Fragment Check"] = bad_results
    df["Overall Match Status"] = bresults
    return df


def filter_excluded_identifications(df):
    """
    Remove identification types that should not be evaluated or exported.

    Excludes:
        - formate-adduct rows
        - all FAHFA rows
        - DG sodium-adduct rows

    Returns:
        filtered_df
        removed_df
    """

    filtered = df.copy()

    filtered["Ontology"] = (
        filtered["Ontology"]
        .astype(str)
        .str.strip()
    )

    filtered["Clean Adduct"] = (
        filtered["Adduct type"]
        .apply(lambda value: fix_adduct(value, 0, 1))
    )

    filtered["Removal Reason"] = ""

    filtered.loc[
        filtered["Clean Adduct"].eq("[M+HCOO]"),
        "Removal Reason"
    ] = "Formate adduct excluded"

    filtered.loc[
        filtered["Ontology"].eq("FAHFA"),
        "Removal Reason"
    ] = "FAHFA excluded"

    filtered.loc[
        (
            filtered["Ontology"].eq("DG")
            & filtered["Clean Adduct"].eq("[M+Na]")
        ),
        "Removal Reason"
    ] = "DG sodium adduct excluded"

    removed_df = filtered[
        filtered["Removal Reason"] != ""
    ].copy()

    filtered_df = filtered[
        filtered["Removal Reason"] == ""
    ].copy()

    print("\nRows removed:")
    print(
        removed_df["Removal Reason"]
        .value_counts()
        .to_string()
    )

    print(
        f"\nRows remaining after exclusions: "
        f"{len(filtered_df)}"
    )

    filtered_df = filtered_df.drop(
        columns=["Clean Adduct", "Removal Reason"]
    )

    removed_df = removed_df.drop(
        columns=["Clean Adduct"]
    )

    return filtered_df, removed_df

def print_check_results(df):
    # For each unique class, print the percentage of Overall Match Status that are True
    for lipid_class in df["Ontology"].unique():
        class_df = df[df["Ontology"] == lipid_class]
        total = len(class_df)
        matched = len(class_df[class_df["Overall Match Status"] == True])
        percentage = (matched / total) * 100 if total > 0 else 0
        print(f"{lipid_class}: {matched}/{total} ({percentage:.2f}%) matched required fragments.")

    print("Overall Rate: ", len(df[df["Overall Match Status"] == True]), "/", len(df), "(",
          (len(df[df["Overall Match Status"] == True]) / len(df)) * 100, "%)")


if __name__ == "__main__":
    filepath = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline\Combined_DDA_Full.xlsx"

    os.chdir(os.path.dirname(filepath))
    df = pd.read_excel(filepath)

    # Remove excluded identification types before fragment checking.
    # This keeps FAHFA and DG [M+Na] out of the validation results,
    # output spreadsheet, and all downstream analyses.
    df, removed_df = filter_excluded_identifications(df)

    df = set_tolcol(df)

    df = set_basic_tail_names(df, "Metabolite name")

    df = check_fragments(
        df,
        plot_individual=False
    )

    # Initial descriptive status
    df["Identification Status"] = df["Overall Match Status"].map(
        {
            True: "Confirmed",
            False: "Failed"
        }
    )

    df = flag_weird_lipids(
        df,
        max_unsaturation_by_class={
            "PC": 6,
            "PE": 6,
            "PG": 6,
            "PI": 6,
            "PS": 6,
            "PA": 6,
            "SM": 4,
            "EtherPC": 6,
            "HexCer": 8
        },
        sm_min_backbone_carbons=16,
        sm_max_backbone_carbons=20,
        sm_max_acyl_carbons=26,
        sm_max_backbone_unsaturation=2,
        sm_max_acyl_unsaturation=3,
        sm_min_total_carbons=30,
        sm_max_total_carbons=46,
        etherpc_min_chain_carbons=24
    )

    # Restore PC [M+H]+ rows that contain the diagnostic m/z 184 ion.
    # These are retained as valid PC-class IDs, but with caution because
    # positive mode alone does not confirm the individual fatty-acyl tails.
    clean_adduct = df["Adduct type"].apply(
        lambda x: fix_adduct(x, 0, 1)
    )

    positive_pc_mask = (
            (df["Ontology"] == "PC")
            & (clean_adduct == "[M+H]")
            & df["Matched Required Fragments"]
            .astype(str)
            .str.contains(
        r"HG\(184\): Yes",
        regex=True,
        na=False
    )
    )

    df.loc[
        positive_pc_mask,
        "Overall Match Status"
    ] = True

    df.loc[
        positive_pc_mask,
        "Identification Status"
    ] = "PC class confirmed in positive mode; retain with caution"

    removed_df.to_excel(
        "Removed_Identifications.xlsx",
        index=False
    )
    # Filter to Ontology == TG
    # df = df[df["Ontology"] == "DG"]

    write_colored_excel(df, "Overall Match Status", "Combined_Fragment_Assignment_Colored.xlsx",
                        color_map={True: "lightgreen", False: "lightcoral"})