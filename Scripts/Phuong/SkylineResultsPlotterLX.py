import os
import PublicScripts.Lipids.SkylineResultsPlotter as srp
import pandas as pd
import io
import matplotlib.pyplot as plt
import numpy as np

israw = r"""
Class,Lipid,Name 2,Concentration (ug/mL),m/z (+),m/z (-),RT
PC,PC 15:0-18:1(+[2]H7),IS PC 15:0-18:1(+[2]H7),1,779.62,777.62,9.2
PG,PG 15:0-18:1(+[2]H7),IS PG 15:0-18:1(+[2]H7),1,768.5,766.5,7.8
"""
equisplash_df = pd.read_csv(io.StringIO(israw))
equisplash_df["Mass"] = (equisplash_df["m/z (+)"] + equisplash_df["m/z (-)"]) / 2.
# Calculate conc in uM from ug/mL using mass
equisplash_df["Concentration (uM)"] = equisplash_df["Concentration (ug/mL)"] / equisplash_df["Mass"]


file = r"Z:\Group Share\Phuong\Molecule Transition Results.csv"
os.chdir(os.path.dirname(file))

df = pd.read_csv(file)

# Drop all where the adduct doesn't contain -
# df = df[df["Precursor Adduct"].str.contains("-")]

# Map precursor adduct of [M58.1989-] to [M57.9855-]
df["Precursor Adduct"] = df["Precursor Adduct"].replace({"[M58.1989-]": "[M57.9855-]"})

normdf, replicates = srp.sum_transitions(df, mode="Products", drop_IS=True, normalize_IS=True,
                                         normalize_TIC=True, customisdf=equisplash_df)


print(normdf)
print(normdf.columns)
# Plot intensity for each row versus replicates for each lipid class
for i, row in normdf.iterrows():
    lipid_class = row["Molecule List Name"]
    name = row["Molecule"]
    intensities = [row[f"{r}"] for r in replicates]
    plt.plot(np.arange(0, len(intensities)), intensities, marker="o", label=f"{name}")
    plt.xticks(np.arange(0, len(intensities)), replicates, rotation=90)
    plt.xlabel("Replicate")
    plt.ylabel("Normalized Intensity")

plt.tight_layout()
plt.show()


normdf.to_excel("SkylineResults_Summed.xlsx", index=False)

