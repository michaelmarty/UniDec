import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib import rcParams

# Font Parameters
#rcParams['font.size']=28
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

# Path to the excel file
path = "Z:\\Group Share\\Julia\\viroporin detergent screen2.xlsx"
# Read in the Combined Sheet
df = pd.read_excel(path)

print(df)

# Parse the data
n = np.array([7, 6, 5, 4, 3, 2, 1])
data = df.iloc[:, 3:3 + len(n)].to_numpy().transpose()

# Normalize the data, NOTE: Comment this out if you don't want to normalize
data = data / np.max(data, axis=0)

# Start to assemble the labels
Detergent = df["Detergent"].to_numpy() + " "
Protein = df["Protein"].to_numpy() + " "

#PC = df["PC"].to_numpy() + " "
#PC = plt.imshow(data, cmap="Blues")



#PG = df["PG"].to_numpy() + " "
#PG = plt.imshow(data, cmap="Greens")
#conc = df["Conc"].to_numpy().astype(str)

# Add the three strings together
label = Protein + Detergent

# Make the plots
plt.figure()
# Generate the squares
plt.imshow(data, cmap="Purples")
# Make the grid
plt.xticks(np.arange(-0.5, len(df) - 0.5, 1), labels=label, rotation=90)
plt.yticks(np.arange(-0.5, len(n) - 0.5, 1), labels=n.astype(str))
plt.grid()

plt.tight_layout()
plt.savefig("Z:\\Group Share\\Julia\\viroporin detergent screen2fixed.pdf")
plt.show()