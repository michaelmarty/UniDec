import os
from Lipidomics_Stats_Tests import *

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_Brain_Classes.xlsx"
path = "NA_Ecoli_Classes.xlsx"

df = pd.read_excel(path, sheet_name="Classes")
tdf = pd.read_excel(path, sheet_name="Tail Lengths")
udf = pd.read_excel(path, sheet_name="Tail Unsaturation")

columns = df.columns[1:]

# Strip off the class from the name
classes = df["Class Name"]

# Get unique classes while preserving order
indexes = np.unique(classes, return_index=True)[1]
uclasses = np.array([classes[index] for index in sorted(indexes)])

output_data = []

for col in columns:
    splits = col.split("_")
    temp = splits[0]
    belt = splits[1]

    print(temp, belt)

    for j, c in enumerate(uclasses):
        intensity = df[col][j]
        tl = tdf[col][j]
        tu = udf[col][j]
        output_data.append([temp, belt, c, intensity, tl, tu])

outdf = pd.DataFrame(output_data,
                     columns=["Temperature", "Belt", "Class Name", "Intensity", "Tail Length", "Tail Unsaturation"])

writer = pd.ExcelWriter(path, engine='openpyxl', mode='a', if_sheet_exists="replace")
outdf.to_excel(writer, index=False, sheet_name="Reformat")
writer.close()
