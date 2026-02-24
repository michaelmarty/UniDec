import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === INPUT FILE ===
file = r"D:\UTAustin Summer 2025\20251026_EC_D1(-)_native_.txt"
label = "D1 (-) Zero-Charge Spectrum"

# === FONT SETTINGS ===
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10
})

# === COLOR ===
color = "#2A788EFF"  # blue-green (color-blind friendly)

# === LOAD DATA ===
df = pd.read_csv(file, delim_whitespace=True, comment="#", header=None)
df.columns = ["Mass_Da", "Intensity"]

# === CLEAN DATA ===
# Remove any rows where intensity is NaN or negative
df = df.dropna()
df = df[df["Intensity"] >= 0]

# Trim leading zeros — keeps only region with signal
nonzero_indices = df.index[df["Intensity"] > 0]
if len(nonzero_indices) > 0:
    start_idx = max(nonzero_indices[0] - 50, 0)  # small buffer before first signal
    df = df.iloc[start_idx:]

# Convert to kDa and normalize to 100%
df["Mass_kDa"] = df["Mass_Da"] / 1000
df["Percent"] = (df["Intensity"] / df["Intensity"].max()) * 100

print(f"✅ Loaded {len(df)} data points. Max intensity = {df['Intensity'].max():.2f}")

# === PLOT ===
fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
ax.plot(df["Mass_kDa"], df["Percent"], color=color, linewidth=1.5, label=label)

# === AXIS STYLE ===
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Focus on the meaningful range automatically
mass_min = df["Mass_kDa"].min()
mass_max = df["Mass_kDa"].max()
ax.set_xlim(mass_min - 0, mass_max + 10)
ax.set_ylim(0, 105)

ax.set_xlabel("Mass (kDa)")
ax.set_ylabel("Relative Intensity (%)")

# === LEGEND ===
ax.legend(frameon=False, loc="upper right")

# === BACKGROUND & LAYOUT ===
fig.patch.set_alpha(0.0)
ax.set_facecolor("none")
plt.tight_layout()

plt.show()
