import pandas as pd
import os
from io import StringIO
from PublicScripts.Lipids.LipidFragPrediction import *
from PublicScripts.Lipids.MatchPlotting import *
from PublicScripts.Lipids.LipidFunctions import parse_spec_string, fix_adduct, match_frag_to_spec, set_tails, spec_to_nl, write_colored_excel
import numpy as np
import re
import ms_entropy as me

# Define CSV

lipid_data = """
Class	Adduct	Required Frags
Cer\t[M+H]	FA1(+O), FA2(+O)
Cer\t[M-H]	FA1(+O), FA2(+O)
Cer\t[M+CH3COO]	FA1(+O), FA2(+O)
CL\t[M-H]	FA1(+O), FA2(+O), FA3(+O), FA4(+O), HG(152)
CoQ	[M+H]	HG(197)
DG	[M+NH4]	FA1(+O), FA2(+O)
FAHFA	[M-H]	FA1(+O)
GM3	[M-H]	HG(290), HGNL(-290)
HexCer	[M+H]	HGNL(162), FA1(+O)
HexCer	[M-H]	HG, FA1(+O)
HexCer	[M+CH3COO]	HG, FA1(+O)
LPA	[M-H]	HG(152)
LPC	[M+H]	HG(184)
LPC	[M+CH3COO]	HG(224), FA1(+O)
LPC O	[M+CH3COO]	HG(168)
LPC O	[M+H]	HG(184)
LPE	[M-H]	HG(196), FA1(+O)
LPG	[M-H]	HG(152), FA1(+O)
LPI	[M-H]	HG(152), FA1(+O), HG(241)
LPS	[M-H]	HG(152), FA1(+O)
PA	[M-H]	HG(152), FA1(+O), FA2(+O)
PC	[M+H]	HG(184)
PC	[M+CH3COO]	HG(224), HGNL(74), FA1(+O), FA2(+O)
EtherPC	[M+CH3COO]	HG(224), FA1(+O), FA2(+O), FA1(+O-C2H4), FA2(+O-C2H4) 
EtherPC	[M+H]	HG(184)
PE	[M+H]	HGNL(141), FA1(+OH-OH), FA2(+OH-OH)
PE	[M-H]	HG(196), FA1(+O), FA2(+O)
EtherPE	[M+H]	HGNL(141), FA1(+OH), FA2(+OH)
EtherPE	[M-H]	FA1(+O), FA2(+O), NLFA1(-H), NLFA2(-H)
EtherPE [M+CH3COO]	FA1(+O), FA2(+O)
PG	[M-H]	HG(152), FA1(+O), FA2(+O)
PI	[M-H]	HG(241), HG(297), FA1(+O), FA2(+O)
PS	[M-H]	HGNL(87), FA1(+O), FA1+GP, FA2+GP
SM	[M+CH3COO]	HGNL(74), HG(168)
SM	[M+H]	HG(184)
TG	[M+NH4]	NLFA1(+O+H2O), NLFA2(+O+H2O), NLFA3(+O+H2O)
"""
# Read the CSV data into a pandas DataFrame

lipid_df = pd.read_csv(StringIO(lipid_data), sep="\t", index_col=False)
# Strip whitespace from classes
lipid_df["Class"] = lipid_df["Class"].str.strip()


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
    # Get row from hgdf_default that matches lipid_class and adduct
    row = hgdf_default[(hgdf_default["CompoundClass"] == lipid_class) & (hgdf_default["Adduct"] == fix_adduct(adduct, 1, 0))]
    if row.empty:
        print(f"No matching row found for {lipid_class} with adduct {adduct}.")
        return None

    if "HGNL" in fragname:
        vals = row["HeadgroupNL_mz"].values[0]
    elif "HG" in fragname:
        vals = row["HeadgroupFragment_mz"].values[0]
    else:
        raise ValueError(f"Fragment name {fragname} not recognized.")

    if vals == "":
        print(f"No values found for {lipid_class} with adduct {adduct} and fragment {fragname}.")
        return None
    vals = vals.split(",")
    fragnum = fragname.split("(")[-1].split(")")[0] if "(" in fragname else None
    if fragnum is None:
        print(f"No fragment number specified in {fragname}.")
        return None

    for v in vals:
        if v.startswith(fragnum):
            return float(v)

    print(f"No matching fragment number {fragnum} found in values {vals} for {lipid_class} with adduct {adduct}.")
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


def check_fragments(df, ldf=None, plot_individual=False):
    if ldf is None:
        ldf = lipid_df

    # Loop through each row in the input DataFrame
    results = []
    bresults = []
    for index, row in df.iterrows():
        lipid_class = row["Ontology"].strip()
        adduct = row["Adduct type"]
        mz = row["Reference m/z"]
        # Clean adduct
        adduct = fix_adduct(adduct, 0, 1)

        # Match ldf to get required fragments
        required_frags = ldf[(ldf["Class"] == lipid_class) & (ldf["Adduct"] == adduct)]["Required Frags"].values

        if len(required_frags) == 0:
            print(f"No required fragments found for {lipid_class} with adduct {adduct}.")
            results.append("No required fragments found for this type")
            bresults.append(False)
            continue

        tol = row["Tolerance"]
        spec = parse_spec_string(row["MS/MS spectrum"], threshold=0.001, norm=True)
        ref_spec = parse_spec_string(row["Ref Spec"], threshold=0.001, norm=True)

        if len(spec) == 0 or len(ref_spec) == 0:
            print(f"Empty spectrum for {row['Metabolite name']}. Skipping.")
            results.append("Empty spectrum")
            bresults.append(False)
            continue

        matches = ""
        mbool = False

        for frags in required_frags:
            matches_sub = ""
            mbool_sub=True
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
                        matched_correctly, match_bool = get_match_result(spec_to_nl(spec, mz), spec_to_nl(ref_spec, mz), frag, tol)
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
                        matched_correctly, match_bool = get_match_result(spec_to_nl(spec, mz), spec_to_nl(ref_spec, mz), frag, tol)
                    else:
                        matched_correctly, match_bool = get_match_result(spec, ref_spec, frag, tol)
                    mbool_sub = mbool_sub and match_bool

                    print(f"FA fragment for {lipid_class} with adduct {adduct}: {frag}", matched_correctly)

                    matches_sub += f"{f}: {matched_correctly}; "

            if mbool_sub:
                mbool = mbool_sub
                matches = matches_sub
                continue

            if matches=="":
                matches+=matches_sub

        # print("Score based on Spectral Entropy:", spectral_entropy_score(spec, ref_spec, tol=tol))
        #
        if plot_individual and not mbool:
            butterfly_plot(spec, ref_spec, title=row["Metabolite name"])

        results.append(matches)
        bresults.append(mbool)


    df["Matched Required Fragments"] = results
    df["Overall Match Status"] = bresults
    return df

def print_check_results(df):
    # For each unique class, print the percentage of Overall Match Status that are True
    for lipid_class in df["Ontology"].unique():
        class_df = df[df["Ontology"] == lipid_class]
        total = len(class_df)
        matched = len(class_df[class_df["Overall Match Status"] == True])
        percentage = (matched / total) * 100 if total > 0 else 0
        print(f"{lipid_class}: {matched}/{total} ({percentage:.2f}%) matched required fragments.")

    print("Overall Rate: ", len(df[df["Overall Match Status"] == True]), "/", len(df), "(", (len(df[df["Overall Match Status"] == True]) / len(df)) * 100, "%)")



if __name__ == "__main__":
    filepath = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline\Lipid Checker\Combined_DDA_Full.xlsx"
    os.chdir(os.path.dirname(filepath))
    df = pd.read_excel(filepath)

    # Filter to Ontology == PC
    df = df[df["Ontology"] == "EtherPE"]
    # print("Filtered PC DataFrame:")

    df = set_tolcol(df)
    df = set_basic_tail_names(df, "Metabolite name")

    # df = check_fragments(df, plot_individual=True)
    df = check_fragments(df, plot_individual=False)
    print_check_results(df)

    # df = assign_df_fragments(df, tolcol="Tolerance", auto_tolcol=True)

    # Write df to test output
    # df.to_excel("Combined_Fragment_Assignment.xlsx", index=False)

    write_colored_excel(df, "Overall Match Status", "Combined_Fragment_Assignment_Colored.xlsx", color_map={True: "lightgreen", False: "lightcoral"})

