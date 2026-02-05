import pandas as pd
import os
import numpy as np

file = r"Y:\Data\Lipidomics\Libraries\LipidCreatorSetup.xlsx"

classdf = pd.read_excel(file, sheet_name="Classes")
taildf = pd.read_excel(file, sheet_name="Tails")

# drop all where count is below 6
taildf = taildf[taildf['Count'] >= 6]
tails = list(taildf['Tails'])
print(len(tails))

lipidlist = []

# Loop over all classes
for index, row in classdf.iterrows():
    cls = row['Class']
    ntails = int(row["NFattyAcids"])

    names = []
    if ntails == 0:
        continue
    elif ntails == 1:
        names = tails
    elif ntails == 2:
        for i in range(len(tails)):
            for j in range(i, len(tails)):
                names.append(f"{tails[i]}_{tails[j]}")
    elif ntails == 3:
        for i in range(len(tails)):
            for j in range(i, len(tails)):
                for k in range(j, len(tails)):
                    names.append(f"{tails[i]}_{tails[j]}_{tails[k]}")
    elif ntails == 4:
        for i in range(len(tails)):
            for j in range(i, len(tails)):
                for k in range(j, len(tails)):
                    for l in range(k, len(tails)):
                        names.append(f"{tails[i]}_{tails[j]}_{tails[k]}_{tails[l]}")
    else:
        print("Too many tails, unsupported")
        raise Exception("Too many tails, unsupported")

    for n in names:
        lipidlist.append(f"{cls} {n}")

    print(f"Class: {cls} Total Lipids: {len(names)}")


print(len(lipidlist))

outfile = r"Y:\Data\Lipidomics\Libraries\LipidCreatorAllLipids.csv"
with open(outfile, 'w') as f:
    for l in lipidlist:
        f.write(f"{l}\n")



