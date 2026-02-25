import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# === SETTINGS ===
files = [
   # r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_BIORAD_secstdtry1 001.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1T0minus_1to80_DMPC 001.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to80DMPC.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to72_DMPC.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1t65_DMPC.xlsx",
    r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to60_DMPC.xlsx",
    #r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_maybeD1sfCherryshort_minustag_1to150DMPC 001.xlsx",
    #r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1sfCherry(minus)LONGLINKER_1to150_DMPC 001.xlsx",
r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_sfGFPD1T0minus_1to150_DMPC 001.xlsx",
#r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1halo1toDMPC 001.xlsx",
]

labels = [
    #"SEC Standards",
    "MSP1D1",
    "MSP1D1-sfGFP 1:80",
    "MSP1D1-sfGFP 1:72",
    "MSP1D1-sfGFP 1:65",
    "MSP1D1-sfGFP 1:60",
    "MSP1D1-sfGFP (no BCA 1:150)",
    #"MSP1D1-sfCherry Short",
    #"MSP1D1-sfCherry Long",
#"MSP1D1-sfGFP",
#"MSP1D1-Halo",
]

if len(files) != len(labels):
    raise ValueError("The number of files must match the number of legend labels.")

# === PLOT CONFIG ===
mpl.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'savefig.transparent': True,
    'axes.facecolor': 'none',
    'figure.facecolor': 'none'
})

colors = [
    "#0F0C08",  # charcoal gray
    "#3CAEA3",  # dark blue
    "#2C7BB6",  # warm brown
    "#22A884FF",  # slate gray
    "#7AD151FF",  # olive gray
    "#A8A29E",  # stone
    "#FDE725FF"   # almost black
]

plt.figure(figsize=(6, 4), facecolor='none')

# === READ AND PLOT EACH FILE ===
for i, file_path in enumerate(files):
    ext = os.path.splitext(file_path)[1].lower()

    try:
        df = pd.read_excel(file_path, header=None, engine="openpyxl")
    except Exception as e:
        print(f"❌ Could not read {file_path}: {e}")
        continue

    # Detect header row
    header_row = None
    for idx in range(min(10, len(df))):
        row = df.iloc[idx].astype(str).str.lower().str.strip().values
        if any("ml" in s for s in row) and any("mau" in s for s in row):
            header_row = idx
            break

    if header_row is not None:
        df = pd.read_excel(file_path, header=header_row, engine="openpyxl")
    else:
        df = pd.read_excel(file_path, header=2, engine="openpyxl")

    # Identify ml and mAU columns
    cols_lower = [c.lower().strip() for c in df.columns]
    ml_col, mau_col = None, None
    for j, cname in enumerate(cols_lower):
        if "ml" in cname:
            ml_col = df.columns[j]
        if "mau" in cname:
            mau_col = df.columns[j]
    if ml_col is None or mau_col is None:
        ml_col, mau_col = df.columns[0], df.columns[1]

    df[ml_col] = pd.to_numeric(df[ml_col], errors="coerce")
    df[mau_col] = pd.to_numeric(df[mau_col], errors="coerce")
    df_clean = df[[ml_col, mau_col]].dropna()
    df_clean.columns = ["ml", "mAU"]
    df_clean = df_clean.sort_values("ml").reset_index(drop=True)

    # === Normalize mAU ===
    max_val = df_clean["mAU"].max()
    if max_val != 0:
        df_clean["mAU"] = df_clean["mAU"] / max_val

    # Plot
    plt.plot(
        df_clean["ml"],
        df_clean["mAU"],
        linewidth=1.3,
        label=labels[i],
        color=colors[i % len(colors)],
        alpha=0.9
    )

# === STYLE ===
ax = plt.gca()

# Axis limits
plt.xlim(left=0, right=26)
plt.ylim(0, 1)

# Axis labels
plt.xlabel("Retention Volume (mL)")
plt.ylabel("Normalized Absorbance")

# Remove top and right spines only (keep bottom & left connected)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Keep spines connected to axes
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

# Ticks and legend
plt.legend(frameon=False)
plt.tight_layout()

# Save transparent
plt.savefig("chromatograms_normalized_clean.png", transparent=True, bbox_inches='tight', dpi=300)
plt.show()

print("✅ Normalized chromatogram (connected axes) saved as 'chromatograms_normalized_clean.png'")