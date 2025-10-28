import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib import rcParams

#rcParams['font.size']=28
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

path = "Z:\\Group Share\\Julia\\viroporin detergent screen.xlsx"

df = pd.read_excel(path)

print(df)

plt.figure()
n = np.array([6, 5, 4, 3, 2, 1])

data = df.iloc[:, 3:3 + len(n)].to_numpy()

detergent = df["Detergent"].to_numpy() + " "
protein = df["Protein"].to_numpy() + " "
#conc = df["Conc"].to_numpy().astype(str)

# Add the three strings together
# pepconc = np.char.add(np.char.add(peptides.astype(str), conc.astype(str)), lipid.astype(str))
label = detergent + protein


plt.figure
plt.imshow(data.transpose())
plt.imshow(data, cmap="Blues")
plt.xticks(np.arange(-0.5, len(df) - 0.5, 1), labels=label, rotation=90)
plt.yticks(np.arange(-0.5, len(n) - 0.5, 1), labels=n.astype(str))
plt.grid()

plt.tight_layout()

plt.show()
