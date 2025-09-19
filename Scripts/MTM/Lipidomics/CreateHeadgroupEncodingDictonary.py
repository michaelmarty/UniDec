import pandas as pd
import os
import numpy as np

libdir = r"Y:\Data\Lipidomics\Libraries"
os.chdir(libdir)

nfile = "neglib.csv"
pfile = "poslib.csv"
ndf = pd.read_csv(os.path.join(libdir, nfile))
pdf = pd.read_csv(os.path.join(libdir, pfile))

print("Neg DF shape:", ndf.shape)
print("Pos DF shape:", pdf.shape)

nheads = np.unique(ndf["compoundclass"])
pheads = np.unique(pdf["compoundclass"])
print("Neg Headgroups:", len(nheads))
print(nheads)
print("Pos Headgroups:", len(pheads))
print(pheads)

# Combine
allheads = set(nheads).union(set(pheads))
print("All Headgroups:", len(allheads))
print(allheads)

# Add Cer to allheads
allheads.add("Cer")

# Sort and create dictionary
sortedheads = sorted(list(allheads))
headgroup_dict = {head: idx for idx, head in enumerate(sortedheads)}
print("Headgroup Dictionary:", headgroup_dict)
# Save to CSV
headgroup_df = pd.DataFrame(list(headgroup_dict.items()), columns=["Headgroup", "Encoding"])
headgroup_df.to_csv("Headgroup_Encoding_Dictionary.csv", index=False)
