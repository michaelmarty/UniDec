import pandas as pd
import os
from io import StringIO
from PublicScripts.Lipids.LipidFragPrediction import *
from PublicScripts.Lipids.MatchPlotting import *
from PublicScripts.Lipids.LipidFunctions import parse_spec_string, fix_adduct, match_frag_to_spec
import numpy as np
import ms_entropy as me

# Define CSV

lipid_data = """
Class	Adduct	Required Frags
Cer	[M+H]	FA1, FA2
Cer	[M-H]	FA1, FA2
Cer	[M+CH3COO]	FA1, FA2
CL 	[M-H]	FA1, FA2, FA3, FA4
CoQ	[M+H]	HG
DG	[M+NH4]	FA1, FA2
FAHFA	[M-H]	FA1
GM3	[M-H]	HG
HexCer	[M+H]	HG, FA1
HexCer	[M-H]	HG, FA1
HexCer	[M+CH3COO]	HG, FA1
LPA	[M-H]	HG
LPC	[M+H]	HG
LPC O	[M+CH3COO]	HG
LPC	[M+CH3COO]	HG, FA1
LPC O	[M+H]	HG
LPE	[M-H]	HG, FA1
LPG	[M-H]	HG, FA1
LPI	[M-H]	HG, FA1
LPS	[M-H]	HG, FA1
PA	[M-H]	HG, FA1, FA2
PC	[M+H]	HG(184)
PC O	[M+CH3COO]	HG(224), FA1, FA2
PC	[M+CH3COO]	HG(224), FA1, FA2
PC O	[M+H]	HG(184)
PE	[M+H]	HG, FA1, FA2
PE	[M-H]	HG, FA1, FA2
PE O	[M+H]	HG, FA1, FA2
PE O	[M-H]	HG, FA1, FA2
PG	[M-H]	HG, FA1, FA2
PI	[M-H]	HG, FA1, FA2
PS	[M-H]	HG, FA1, FA2
SM	[M+CH3COO]	HG
SM	[M+H]	HG(184)
TG	[M+NH4]	FA1, FA2, FA3

"""
# Read the CSV data into a pandas DataFrame

lipid_df = pd.read_csv(StringIO(lipid_data), sep="\t")

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



def check_fragments(df, ldf=None):
    if ldf is None:
        ldf = lipid_df

    # Loop through each row in the input DataFrame
    results = []
    for index, row in df.iterrows():
        lipid_class = row["Ontology"]
        adduct = row["Adduct type"]
        # Clean adduct
        adduct = fix_adduct(adduct, 0, 1)

        matched_correctly = "Unknown"

        # Match ldf to get required fragments
        required_frags = ldf[(ldf["Class"] == lipid_class) & (ldf["Adduct"] == adduct)]["Required Frags"].values
        if len(required_frags) == 0:
            print(f"No required fragments found for {lipid_class} with adduct {adduct}.")
            continue


        tol = row["Tolerance"]
        spec = parse_spec_string(row["MS/MS spectrum"], threshold=0.01, norm=True)
        ref_spec = parse_spec_string(row["Ref Spec"], threshold=0.01, norm=True)

        for f in required_frags[0].split(","):
            f = f.strip()
            if "HG" in f:
                frag = get_hg_frag(lipid_class, adduct, f)
                matched_frag = match_frag_to_spec(spec, frag, tol)
                matched_ref = match_frag_to_spec(ref_spec, frag, tol)

                if matched_frag and matched_ref:
                    matched_correctly = "Yes"
                elif not matched_frag and matched_ref:
                    matched_correctly = "No"
                elif matched_frag and not matched_ref:
                    matched_correctly = "Yes but Ref Failed"
                else:
                    matched_correctly = "No and Ref Failed"

                print(f"HG fragment for {lipid_class} with adduct {adduct}: {frag}", matched_frag, matched_ref)

        # print("Score based on Spectral Entropy:", spectral_entropy_score(spec, ref_spec, tol=tol))
        #
        # butterfly_plot(spec, ref_spec, title=row["Metabolite name"])

        results.append(matched_correctly)

    df["Matched Required Fragments"] = results
    return df





if __name__ == "__main__":
    filepath = r"Z:\Group Share\Annika\DDA_Compare\mzmine\Combined_DDA_Full.xlsx"
    os.chdir(os.path.dirname(filepath))
    df = pd.read_excel(filepath)

    # Filter to Ontology == PC
    df = df[df["Ontology"] == "PC"]
    print("Filtered PC DataFrame:")

    df = set_tolcol(df)

    df = check_fragments(df)

    # df = assign_df_fragments(df, tolcol="Tolerance", auto_tolcol=True)

    # Write df to test output
    df.to_excel("Combined_Fragment_Assignment.xlsx", index=False)
