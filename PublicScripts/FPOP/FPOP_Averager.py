import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl

# Some Plotting Parameters
mpl.use("WxAgg")
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 11
rcParams['font.sans-serif'] = "Arial"

# Set the file and directory
topdir = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\For Public Repository\\FPOP Data"
os.chdir(topdir)
file = "RawExtractedPercents.xlsx"

# Read File
df = pd.read_excel(file)

# Set List of Lipids and Ratios to Consider
lipids = ["DLPC", "DMPG", "DMPC", "POPC", "DPPC", "Ecoli", "Free"]
ratios = [3, 9, 18]

# Set Colors
colors = ["pink", "b", "k", "orange", "purple", "g", "red"]

# Loop over all lipids and ratios to get averages and std devs
output = []
for i, l in enumerate(lipids):
    for j, r in enumerate(ratios):
        # Get rows that match lipid type and ratios
        df2 = df[df["Lipid Type"] == l]
        df2 = df2[df2["Ratio"] == r]

        # Exctract relative ox value
        ox = df2["+Ox"].to_numpy()
        # Calculate, n, mean, std dev
        n = len(ox)
        if n > 0:
            m = np.mean(ox)
            s = np.std(ox, ddof=1)
            out = ["Da", l, r, m, s, n]
        else:
            out = ["Da", l, r, 0, 0, n]
        # Add it to the list
        output.append(out)

        # Create Bar Plot
        plt.bar(i * 3 + j, m, yerr=s, color=colors[i])

# Write Outputs
df3 = pd.DataFrame(output)
df3.columns = ["Peptide", "Lipid", "Ratio", "Mean", "SD", "N"]
print(df3)
df3.to_excel("Averager_Output.xlsx")

# Finish Plotting
plt.ylim(0, 75)
plt.show()
