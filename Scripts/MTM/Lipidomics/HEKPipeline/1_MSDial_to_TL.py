import os
from PublicScripts.Lipids.LipidFunctions import *

# Set the working directory
topdir = r"Z:\Group Share\Shivam\Reducing_FP\Finished results"
os.chdir(topdir)

# Set the results file names, should be the peak_master files from MS-DIAL
negfile = "Stellar_neg.csv"
posfile = "Stellar_pos.csv"

# Set the XML files with reference spectra
negxml = "Stellar_neg.xml"
posxml = "Stellar_pos.xml"

# Load the data
negdf = pd.read_csv(negfile)
posdf = pd.read_csv(posfile)

# Cleanup the data
# Drop if comment has "Low Quality" in it (case-insensitive). If no 'Comment' column exists, leave unchanged.
if 'Comment' in posdf.columns:
    posdf = posdf[~posdf['Comment'].str.contains('Low Quality', case=False, na=False)]
    negdf = negdf[~negdf['Comment'].str.contains('Low Quality', case=False, na=False)]

# Drop CL in positive mode
posdf = posdf[~posdf["Ontology"].str.startswith("CL")]

# Sort by Ontology and Name
posdf = posdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)
negdf = negdf.sort_values(by=["Ontology", "Metabolite name"]).reset_index(drop=True)

# Drop if simple dot product is below 0.45
# posdf = posdf[posdf["Simple dot product"] > 0.45]
# negdf = negdf[negdf["Simple dot product"] > 0.45]
#
# make_plots(posdf, negdf)

# Load the reference spectra from XML
refnegdf = parse_xml(negxml)
refposdf = parse_xml(posxml)

# Add reference spectra to the dataframes
negdf = add_ref_from_xml(negdf, refnegdf)
posdf = add_ref_from_xml(posdf, refposdf)

# Merge positive and negative dataframes
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
# postlfileout = posfile.replace(".csv", "_TransitionList.csv")
# postldf.to_csv(postlfileout, index=False)

negtldf = msdial_to_transitionlist(negdf, cutoff=cutoff, spectrumkey=spectrumkey)
# negtlfileout = negfile.replace(".csv", "_TransitionList.csv")
# negtldf.to_csv(negtlfileout, index=False)

# Merge pos and neg Transition Lists
tldf = pd.concat([postldf, negtldf])
# Write to csv
tldf.to_csv("Combined_TransitionList.csv", index=False)
# Print number of unique Molecule Names
print("Number of unique Molecule Names:", tldf["Molecule Name"].nunique())