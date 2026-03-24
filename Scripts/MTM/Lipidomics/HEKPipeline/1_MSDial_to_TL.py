import os
from PublicScripts.Lipids.LipidFunctions import *

# Set the working directory
topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults"
os.chdir(topdir)

# Set the results file names, should be the peak_master files from MS-DIAL
negfile = "Stellar_neg.csv"
posfile = "Stellar_pos.csv"

# Set the XML files with reference spectra
negxml = "Stellar_neg.xml"
posxml = "Stellar_pos.xml"

mergeddf, posdf, negdf = cleanup_and_merge_msdial(posfile, negfile, posxml, negxml, rttol=0.2)

print("posdf columns:", posdf.columns.tolist())
print("negdf columns:", negdf.columns.tolist())

print("'Ref Spec' in posdf:", "Ref Spec" in posdf.columns)
print("'Ref Spec' in negdf:", "Ref Spec" in negdf.columns)
# Cleanup the data
# Drop if comment has "Low Quality" in it (case-insensitive). If no 'Comment' column exists, leave unchanged.
if 'Comment' in posdf.columns:
    posdf = posdf[~posdf['Comment'].str.contains('Low Quality', case=False, na=False)]
    negdf = negdf[~negdf['Comment'].str.contains('Low Quality', case=False, na=False)]

# Drop CL and PS in positive mode
posdf = posdf[~posdf["Ontology"].str.startswith("CL")]
posdf = posdf[~posdf["Ontology"].str.startswith("PS")]
print("POS adducts:")
print(posdf["Adduct type"].value_counts())

print("\nNEG adducts:")
print(negdf["Adduct type"].value_counts())
# Drop if simple dot product is below 0.45
# posdf = posdf[posdf["Simple dot product"] > 0.45]
# negdf = negdf[negdf["Simple dot product"] > 0.45]
#
# make_plots(posdf, negdf)

# Save to output
# posfileout = posfile.replace(".csv", "_Ref.csv")
# posdf.to_csv(posfileout, index=False)
# negfileout = negfile.replace(".csv", "_Ref.csv")
# negdf.to_csv(negfileout, index=False)

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

# Correct Explicit Retention Time by adding 0.5 minutes
if "Explicit Retention Time" in tldf.columns:
    tldf["Explicit Retention Time"] = tldf["Explicit Retention Time"] + 0.5

# Write to csv
tldf.to_csv("Combined_TransitionList.csv", index=False)
# Print number of unique Molecule Names
print("Number of unique Molecule Names:", tldf["Molecule Name"].nunique())

print("postldf rows:", len(postldf))
print("negtldf rows:", len(negtldf))