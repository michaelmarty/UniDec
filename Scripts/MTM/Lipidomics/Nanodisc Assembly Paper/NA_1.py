import os

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_230501_Brain_Neg_TL_Results_cleaned.xlsx"
path2 = "NA_230501_Brain_POS_TL_Results_cleaned.xlsx"


# Load the data
posdf = pd.read_excel(path2, skiprows=1)
negdf = pd.read_excel(path, skiprows=1)

# Rename column molecule name to name
posdf = posdf.rename(columns={"Molecule Name": "Name"})
negdf = negdf.rename(columns={"Molecule Name": "Name"})

namedf = pd.read_excel("NA_230501_Brain_Both_fdata.xlsx")

for i, row in namedf.iterrows():
    sampleid = row["Sample ID"]
    samplename = row["Sample Name"]
    mode = row["mode"]

    if mode == "P":
        posdf = posdf.rename(columns={sampleid: samplename})

    if mode == "N":
        negdf = negdf.rename(columns={sampleid: samplename})



# Add positive and negative column to each data frame
posdf["Pos"] = 1
negdf["Pos"] = 0
posdf["Neg"] = 0
negdf["Neg"] = 1

# Concatenate the data frames
df = pd.concat([posdf, negdf], axis=0)
#df = posdf

df = set_tails(df, "Name")

df = set_class_name(df, "Name")

df = calc_tail_diff(df)

# Loop over column names and change them to lower case
for col in df.columns:
    if "AChR" in col:
        # Convert AChR to AchR
        df = df.rename(columns={col: col.replace("AChR", "AchR")})

#Save the results
df.to_excel("NA_Brain_Both_output.xlsx", index=False)
