import os

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_230502_Ecoli_Neg_TL_Results_cleaned.xlsx"


# Load the data
negdf = pd.read_excel(path, sheet_name="RawData")

# Rename column molecule name to name
negdf = negdf.rename(columns={"Molecule Name": "Name"})

negdf = set_tails(negdf, "Name")

negdf = set_class_name(negdf, "Name")

negdf = calc_tail_diff(negdf)

#Save the results
negdf.to_excel("NA_Ecoli_output.xlsx", index=False)
