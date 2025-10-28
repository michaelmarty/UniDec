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
path = "Z:\\Group Share\\Julia\\viroporin detergent screen.xlsx"

# Read in the Combined Sheet
df = pd.read_excel(path)

print(df)

# Parse the data
n = np.array([0, 1, 2, 3, 4, 5, 6])
ylabel = "Stoichometry"
data = df.iloc[:, 3:3 + len(n)].to_numpy().transpose()

# Normalize the data, NOTE: Comment this out if you don't want to normalize
data = data / np.max(data, axis=0)
data = df.iloc[:, 3:3 + len(n)].to_numpy()
# Start to assemble the labels

Detegent = df["Detergent"].to_numpy() + " "
#peptide = df["peptide"].to_numpy() + " "


# Add the three strings together
#label = Detergent




# Make the plots
plt.figure()
# Generate the squares
plt.imshow(data, cmap="Greens")
# Make the grid
plt.xticks(np.arange(-0.5, len(df) - 0.5, 1), rotation=90)
plt.yticks(np.arange(-0.5, len(n) - 0.5, 1), labels=n.astype(str))
plt.grid()

plt.tight_layout()

plt.show()
