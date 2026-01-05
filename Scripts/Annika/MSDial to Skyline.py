import os
from LipidFunctions import *

# Set the working directory
topdir = r"Z:\Group Share\Annika\IQ-X\UntDDA Tests HEK\New Libs Tests\Manually Filtered\Trans List"
os.chdir(topdir)

# Set the results file names, should be the peak_master files from MS-DIAL
negfile = "PeakMaster_1_2025_12_12_12_57_16.csv"
posfile = "HEKpos_ITIT_Peakmaster.csv"

# Set the XML files with reference spectra
negxml = "SpectrumAndReference_20251211_HEKneg_OTIT_2025_12_12_12_57_19.xml"
posxml = "SpectrumAndReference_20251211_HEKpos_ITIT_2025_12_11_13_48_16.xml"

# Load the data
negdf = pd.read_csv(negfile)
posdf = pd.read_csv(posfile)

make_plots(posdf, negdf)

# Load the reference spectra from XML
refnegdf = parse_xml(negxml)
refposdf = parse_xml(posxml)

# Add reference spectra to the dataframes
negdf = add_ref_from_xml(negdf, refnegdf)
posdf = add_ref_from_xml(posdf, refposdf)

# # Merge positive and negative dataframes to ensure no duplicates
rttol = 0.1
posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

# Save to output
posfileout = posfile.replace(".csv", "_Ref.csv")
posdf.to_csv(posfileout, index=False)
negfileout = negfile.replace(".csv", "_Ref.csv")
negdf.to_csv(negfileout, index=False)

# Convert these to TransitionLists
spectrumkey = "Ref Spec" # or "MS/MS spectrum" if experimental spectra are used, if so, a cutoff is recommended
cutoff=0

postldf = msdial_to_transitionlist(posdf, cutoff=cutoff, spectrumkey=spectrumkey)
postlfileout = posfile.replace(".csv", "_TransitionList.csv")
postldf.to_csv(postlfileout, index=False)

negtldf = msdial_to_transitionlist(negdf, cutoff=cutoff, spectrumkey=spectrumkey)
negtlfileout = negfile.replace(".csv", "_TransitionList.csv")
negtldf.to_csv(negtlfileout, index=False)

# Merge pos and neg Transition Lists
tldf = pd.concat([negtldf,postldf])
# Write to csv
tldf.to_csv("Combined_TransitionList_HEK.csv", index=False)
# Print number of unique Molecule Names
print("Number of unique Molecule Names:", tldf["Molecule Name"].nunique())

# # if False:
#     Convert to LipidList format
#     posldf = msdial_to_lipidlist(posdf)
#     posllfileout = posfile.replace(".csv", "_LipidList.csv")
#     posldf.to_csv(posllfileout, index=False)
#
#     negldf = msdial_to_lipidlist(negdf)
#     negllfileout = negfile.replace(".csv", "_LipidList.csv")
#     negldf.to_csv(negllfileout, index=False)
#
#  # Merge pos and neg Lipid Lists
#     lldf = pd.concat([posldf, negldf])
#     lldf = lldf.drop_duplicates()
#
#     lldf = lldf.sort_values(by=["Molecule List Name", "Molecule Name"])
#     lldf.to_csv("Combined_Lipids.csv", index=False)
#     lldf.to_csv("Combined_LipidList.csv", columns=["Full Name Adduct"], index=False, header=False)
#     # Print number of unique Molecule Names
#     print("Number of unique Molecule Names:", lldf["Molecule Name"].nunique())
#
# # Then Run LCexe.py
# # Then Run BadLipids.py to filter out bad lipids