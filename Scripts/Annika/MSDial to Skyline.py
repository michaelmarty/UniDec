import os
from PublicScripts.Lipids.LipidFunctions import *
import pandas as pd

# Set the working directory
topdir = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Full Dataset From SS\Peak Masters & xml"
os.chdir(topdir)

# Set the results file names, should be the peak_master files from MS-DIAL
# negfile = "Negative_updated_PeakMaster_1_2026_01_23_14_07_18.csv"
posfile = "Positive_updated_PeakMaster_7_2026_01_23_14_30_21.csv"

# Set the XML files with reference spectra
# negxml = "Negative_updated_SpectrumAndReference_AlignmentResult_2026_01_21_15_54_37_2026_01_23_14_07_18.xml"
posxml = "Positive_updated_SpectrumAndReference_AlignmentResult_2026_01_14_17_58_45_2026_01_23_14_30_22.xml"

# Load the data
# negdf = pd.read_csv(negfile)
posdf = pd.read_csv(posfile)

# make_plots(posdf, negdf)  # optional plotting

# Load the reference spectra from XML
# refnegdf = parse_xml(negxml)
refposdf = parse_xml(posxml)

# Add reference spectra to the dataframes
# negdf = add_ref_from_xml(negdf, refnegdf)
posdf = add_ref_from_xml(posdf, refposdf)
#
# # Merge positive and negative dataframes to ensure no duplicates
# rttol = 0.1
# posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

# Save to output
posfileout = posfile.replace(".csv", "_Ref.csv")
posdf.to_csv(posfileout, index=False)
# negfileout = negfile.replace(".csv", "_Ref.csv")
# negdf.to_csv(negfileout, index=False)

# Convert these to TransitionLists
spectrumkey = "Ref Spec"  # or "MS/MS spectrum" if experimental spectra are used
cutoff = 0

postldf = msdial_to_transitionlist(posdf, cutoff=cutoff, spectrumkey=spectrumkey)
# negtldf = msdial_to_transitionlist(negdf, cutoff=cutoff, spectrumkey=spectrumkey)

# ---- CLEAN SKYLINE ADDUCTS ----
def clean_precursor_adduct(df, col="Precursor Adduct"):
    """
    Cleans the Precursor Adduct column in a transition list so Skyline
    does not invent weird adducts like [M.00022+H]+.
    """
    if col not in df.columns:
        raise ValueError(f"No {col} column found. Columns: {df.columns.tolist()}")

    df[col] = (
        df[col].astype(str)
        .str.replace(r"\.\d+", "", regex=True)  # remove tiny decimal offsets
        .str.strip()
    )
    return df

# Apply to positive and negative transition lists
postldf = clean_precursor_adduct(postldf)
# negtldf = clean_precursor_adduct(negtldf)

# Save individual transition lists
postlfileout = posfile.replace(".csv", "_TransitionList.csv")
postldf.to_csv(postlfileout, index=False)
#
# negtlfileout = negfile.replace(".csv", "_TransitionList.csv")
# negtldf.to_csv(negtlfileout, index=False)

# # Merge pos and neg Transition Lists
# tldf = pd.concat([postldf, negtldf])

# Optional: drop duplicates if any
# tldf = tldf.drop_duplicates()

# Write combined CSV
tldf.to_csv("HEK_Pos_TransitionList.csv", index=False)

# Print number of unique Molecule Names
print("Number of unique Molecule Names:", tldf["Molecule Name"].nunique())

if False:
    # Convert to LipidList format (optional)
    posldf = msdial_to_lipidlist(posdf)
    posllfileout = posfile.replace(".csv", "_LipidList.csv")
    posldf.to_csv(posllfileout, index=False)

    negldf = msdial_to_lipidlist(negdf)
    negllfileout = negfile.replace(".csv", "_LipidList.csv")
    negldf.to_csv(negllfileout, index=False)

    # Merge pos and neg Lipid Lists
    lldf = pd.concat([posldf, negldf])
    lldf = lldf.drop_duplicates()

    lldf = lldf.sort_values(by=["Molecule List Name", "Molecule Name"])
    lldf.to_csv("Combined_Lipids.csv", index=False)
    lldf.to_csv("Combined_LipidList.csv", columns=["Full Name Adduct"], index=False, header=False)

    # Print number of unique Molecule Names
    print("Number of unique Molecule Names:", lldf["Molecule Name"].nunique())

# Then Run LCexe.py
# Then Run BadLipids.py to filter out bad lipids
