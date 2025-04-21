import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

# matplotlib.use('WxAgg')

libcsvpath = "C:\\Data\\Lipidomics\\neglib.csv"
os.chdir(os.path.dirname(libcsvpath))

df = pd.read_csv(libcsvpath)
print(df.keys())

df["compoundclass"].astype("category")
df["compoundcodes"] = df["compoundclass"].astype("category").cat.codes

mz = df["precursor_mz"].to_numpy()
rt = df["retention_time"].to_numpy()
ccs = df["ccs"].to_numpy()

print(df["precursortype"])

exit()

for i, m in enumerate(mz):
    diff = mz - m
    b1 = diff > 0.00002
    b2 = diff < 0.01
    b3 = b1 * b2
    if np.any(b3):
        row = df.iloc[i]
        matchindex = np.argmax(b3)
        matchrow = df.iloc[matchindex]
        if row.formula != matchrow.formula and np.abs(row.retention_time - matchrow.retention_time) < 0.5 and np.abs(row.ccs - matchrow.ccs) < 10:
            print("DIFF", row.precursor_mz - matchrow.precursor_mz, "M1 = ", row.compound_name, row.precursor_mz, "M2 = ", matchrow.compound_name, matchrow.precursor_mz)
            print(row.formula, matchrow.formula, row.precursortype, matchrow.precursortype, matchrow.smiles)


exit()

plt.figure(figsize=(8, 8))
plt.subplot(311)
plt.scatter(mz, rt, c=df.compoundcodes, s=0.5, label=df.compoundclass)
plt.xlabel('m/z')
plt.ylabel('RT')
plt.subplot(312)
plt.scatter(mz, ccs, c=df.compoundcodes, s=0.5, label=df.compoundclass)
plt.xlabel('m/z')
plt.ylabel('CCS')
plt.subplot(313)
plt.scatter(rt, ccs, c=df.compoundcodes, s=0.5, label=df.compoundclass)
plt.xlabel('RT')
plt.ylabel('CCS')
plt.tight_layout()
plt.show()
