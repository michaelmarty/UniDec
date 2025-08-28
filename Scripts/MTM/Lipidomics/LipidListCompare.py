import pandas as pd
import os
from LipidFunctions import remove_specific_problems

topdir = r"Z:\Group Share\Annika\Stellar"
file2 = "Sciex_Combined_LipidList.csv"
file1 = "Thermo Merged LipidList.csv"
os.chdir(topdir)

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

list1 = set(df1["Molecule Name"])
list2 = set(df2["Molecule Name"])
print("File 1 length:", len(list1))
print("File 2 length:", len(list2))
only1 = list1 - list2
only2 = list2 - list1
print("Only in file 1:", len(only1))
print("Only in file 2:", len(only2))
print("Common in both:", len(list1.intersection(list2)))
# Create a merged DF with a column indicating the source file or both
merged = pd.DataFrame(columns=["Molecule List Name", "Molecule Name", "Full Name Adduct", "Source"])
for name in list1.union(list2):
    if name in list1 and name in list2:
        source = "Both"
        row = df1[df1["Molecule Name"] == name].iloc[0]
    elif name in list1:
        source = "Thermo"
        row = df1[df1["Molecule Name"] == name].iloc[0]
    else:
        source = "Sciex"
        row = df2[df2["Molecule Name"] == name].iloc[0]
    merged = pd.concat([merged, pd.DataFrame({
        "Molecule List Name": [row["Molecule List Name"]],
        "Molecule Name": [row["Molecule Name"]],
        "Full Name Adduct": [row["Full Name Adduct"]],
        "Source": [source]
    })], ignore_index=True)


# Sort
merged = merged.sort_values(by=["Molecule List Name", "Molecule Name"])

# Remove some specific problems
merged = remove_specific_problems(merged)

print(len(merged))

merged.to_csv("LipidList_Comparison.csv", index=False)

# Save just the Full Name Adduct to CSV with no header
with open("LipidList_Comparison_Import.csv", 'w') as f:
    for name in merged["Full Name Adduct"]:
        # Replace _ in name with -
        name = name.replace("_", "-")
        f.write(f"{name}\n")