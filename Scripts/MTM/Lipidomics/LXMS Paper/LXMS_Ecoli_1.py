import os

os.chdir("C:\Data\Lipidomics\\AmtBAqpZ")

# Set Paths
path = "230427_Ecoli_Export_Raw.xlsx"

# Load the data
df = pd.read_excel(path, sheet_name="RawData")

# Rename column molecule name to name
df = df.rename(columns={"Molecule Name": "Name"})

namedf = pd.read_csv("230413_Ecoli_LXMS_all_fdata_short.csv")

# Concatenate the data frames

df = set_tails(df, "Name")

df = set_class_name(df, "Name")

df = calc_tail_diff(df)

for i, row in namedf.iterrows():
    sampleid = row["SampleName"]
    samplename = row["Combined"]
    number = row["Sample Number"]
    samplename = samplename + "_" + str(number)
    print(sampleid, samplename)
    df = df.rename(columns={sampleid: samplename})

print(df.columns)
#Save the results
df.to_excel("LXMS_Ecoli_Output.xlsx", index=False)
