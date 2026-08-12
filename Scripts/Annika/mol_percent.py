import pandas as pd
import matplotlib.pyplot as plt

# Read your CSV file
df = pd.read_csv(r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\BB vs DC\Cholesterol_mol%_BB.csv")

# Make sure these columns are numbers
df["Peak Area"] = pd.to_numeric(df["Peak Area"])
df["IS Area"] = pd.to_numeric(df["IS Area"])
df["IS umol"] = pd.to_numeric(df["IS umol"])

# Correct each lipid by its internal standard
df["Corrected umol"] = (df["Peak Area"] / df["IS Area"]) * df["IS umol"]

# Keep only Cholesterol and POPC
df = df[df["Lipid"].isin(["Cholesterol", "POPC"])]

# Reshape so each sample has one Cholesterol value and one POPC value
wide = df.pivot(index="Sample", columns="Lipid", values="Corrected umol")

# Calculate cholesterol mol%
wide["Cholesterol mol%"] = (
    wide["Cholesterol"] / (wide["Cholesterol"] + wide["POPC"])
) * 100

# Print results
print(wide)

plt.figure(figsize=(6, 4))

plt.bar(wide.index, wide["Cholesterol mol%"])

plt.ylabel("Cholesterol mol%")
plt.xlabel("Sample")
plt.title("Cholesterol mol% compared to POPC")

plt.tight_layout()
plt.show()