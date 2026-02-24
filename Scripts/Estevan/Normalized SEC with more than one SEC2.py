import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# === SETTINGS ===
mau_file = r"D:\UTAustin Summer 2025\xlsx files\.5percent_DDM_rep2_mAU.xlsx"
mv_file = r"D:\UTAustin Summer 2025\xlsx files\.5percent_DDM_rep2_mV.xlsx"

mau_label = "AqpZ (mAU)"
mv_label = "AqpZ (mV)"

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
    "#3CAEA3",  # teal
]

plt.figure(figsize=(6, 4), facecolor='none')
ax1 = plt.gca()


# ===== Function to load XLSX and detect columns =====
def load_chrom_file(path, signal_keywords):
    try:
        df = pd.read_excel(path, header=None, engine="openpyxl")
    except Exception as e:
        print(f"❌ Could not read {path}: {e}")
        return None

    # Try to detect header row
    header_row = None
    for idx in range(min(10, len(df))):
        row = df.iloc[idx].astype(str).str.lower().str.strip().values
        if any("ml" in s for s in row) and any(any(key in s for key in signal_keywords) for s in row):
            header_row = idx
            break

    if header_row is not None:
        df = pd.read_excel(path, header=header_row, engine="openpyxl")
    else:
        df = pd.read_excel(path, header=2, engine="openpyxl")

    # Identify columns
    cols_lower = [c.lower().strip() for c in df.columns]
    ml_col, sig_col = None, None

    for j, cname in enumerate(cols_lower):
        if "ml" in cname:
            ml_col = df.columns[j]
        for key in signal_keywords:
            if key in cname:
                sig_col = df.columns[j]

    # fallback → first two columns
    if ml_col is None or sig_col is None:
        ml_col, sig_col = df.columns[0], df.columns[1]

    df[ml_col] = pd.to_numeric(df[ml_col], errors="coerce")
    df[sig_col] = pd.to_numeric(df[sig_col], errors="coerce")

    clean = df[[ml_col, sig_col]].dropna()
    clean.columns = ["ml", "signal"]
    clean = clean.sort_values("ml").reset_index(drop=True)
    return clean


# ===== Load mAU trace =====
df_mau = load_chrom_file(mau_file, ["mau"])
if df_mau is not None:
    df_mau["signal"] /= df_mau["signal"].max()  # normalize
    ax1.plot(df_mau["ml"], df_mau["signal"], color=colors[0],
             label=mau_label, linewidth=1.3)

# ===== Load mV trace (second file, second axis) =====
df_mv = load_chrom_file(mv_file, ["mv"])
ax2 = ax1.twinx()

if df_mv is not None:
    df_mv["signal"] /= df_mv["signal"].max()  # normalize
    ax2.plot(df_mv["ml"], df_mv["signal"], color=colors[1],
             label=mv_label, linewidth=1.3, alpha=0.8)

# === STYLE ===
ax1.set_xlim(0, 26)
ax1.set_ylim(0, 1)
ax2.set_ylim(0, 1)

ax1.set_xlabel("Retention Volume (mL)")
ax1.set_ylabel("Normalized Absorbance (mAU)")
ax2.set_ylabel("Normalized Voltage (mV)")

# Clean spines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_position(('data', 0))
ax1.spines['bottom'].set_position(('data', 0))

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
plt.legend(lines1 + lines2, labels1 + labels2, frameon=False, loc="upper right")

plt.tight_layout()
plt.savefig("chromatograms_mAU_mV.png", transparent=True, dpi=300)
plt.show()

print("✅ Combined chromatogram saved as 'chromatograms_mAU_mV.png'")
