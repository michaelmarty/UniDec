import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

data = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\PosIDs3.csv"
os.chdir("Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX")

df = pd.read_csv(data)

libfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-Corrected2.msp"

# Read library file
with open(libfile, 'r') as f:
    lines = f.readlines()

# Loop through df
rts = []
for index, row in df.iterrows():
    name = row['Metabolite name']
    rt = row['Average Rt(min)']
    refrt = row['Reference RT']
    librt = 0
    try:
        n2 = name.split("|")[1]
    except:
        n2 = name
    try:
        # Replace O2 and O3 with 2O and 3O
        n2 = n2.replace("O2", "2O").replace("O3", "3O")
    except:
        pass
    print(n2)
    # Find name in library
    for i, line in enumerate(lines):
        if line.startswith("NAME: "+n2):
            print(line)
            # Look for RT line
            for j in range(i, min(i+20, len(lines))):
                if lines[j].startswith("RETENTIONTIME:"):
                    librt = float(lines[j].split(":")[1].strip())
                    print("Found RT:", librt)
                    break
            break
    rts.append([name, rt, refrt, librt])

# Convert to DataFrame
rtdf = pd.DataFrame(rts, columns=['Name', 'Average Rt(min)', 'Reference RT', 'Library RT'])
print(rtdf)

# Save to CSV
rtdf.to_csv("RTComparisonPos3Corr.csv", index=False)

exit()


data = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\NegIDs.csv"
os.chdir("Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX")

df = pd.read_csv(data)

libfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp"

# Read library file
with open(libfile, 'r') as f:
    lines = f.readlines()

# Loop through df
rts = []
for index, row in df.iterrows():
    name = row['Metabolite name']
    rt = row['Average Rt(min)']
    refrt = row['Reference RT']
    librt = 0
    try:
        n2 = name.split("|")[1]
    except:
        n2 = name
    try:
        # Replace O2 and O3 with 2O and 3O
        n2 = n2.replace("O2", "2O").replace("O3", "3O")
    except:
        pass
    print(n2)
    # Find name in library
    for i, line in enumerate(lines):
        if line.startswith("NAME: "+n2):
            print(line)
            # Look for RT line
            for j in range(i, min(i+20, len(lines))):
                if lines[j].startswith("RETENTIONTIME:"):
                    librt = float(lines[j].split(":")[1].strip())
                    print("Found RT:", librt)
                    break
            break
    rts.append([name, rt, refrt, librt])

# Convert to DataFrame
rtdf = pd.DataFrame(rts, columns=['Name', 'Average Rt(min)', 'Reference RT', 'Library RT'])
print(rtdf)

# Save to CSV
rtdf.to_csv("RTComparison.csv", index=False)