import os

os.chdir("C:\Data\Lipidomics\\AchR")

# Set Paths
path = "230428_Neg_AchR_IS_fix.csv"
path2 = "230428_Pos_AchR_IS_fix_1.csv"


# Load the data
posdf = pd.read_csv(path2)
negdf = pd.read_csv(path)

# Rename column molecule name to name
posdf = posdf.rename(columns={"Molecule Name": "Name"})
negdf = negdf.rename(columns={"Molecule Name": "Name"})

namedf = pd.read_csv("230403_Brain_AchR_fdata.csv")


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
        df = df.rename(columns={col: col.replace(" ", "")})

for i, row in namedf.iterrows():
    sampleid = row["SampleName"]
    samplename = row["Combined"]
    number = row["Sample Number"]
    samplename = samplename + "_" + str(number)
    print(sampleid, samplename)
    df = df.rename(columns={sampleid: samplename})

print(df.columns)
#Save the results
df.to_excel("LXMS_AchR_Output.xlsx", index=False)
