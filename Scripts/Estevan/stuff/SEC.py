import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# === SETTINGS ===
file_path = r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1T0minus_1to80_DMPC 001.xlsx"

# === PLOT CONFIG ===
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['savefig.transparent'] = True

# === READ FILE SAFELY ===
ext = os.path.splitext(file_path)[1].lower()

if ext not in (".xlsx", ".xls"):
    raise ValueError("Please provide an Excel file (.xlsx or .xls)")

try:
    raw = pd.read_excel(file_path, header=None, engine="openpyxl")
except Exception as e:
    raise RuntimeError(f"Failed to read Excel file: {e}")

# === FIND THE HEADER ROW ===
header_row = None
for idx in range(min(10, len(raw))):  # look through first 10 rows
    row = raw.iloc[idx].astype(str).str.lower().str.strip().values
    if any("ml" == s or s.startswith("ml") for s in row) and any("mau" == s or s.startswith("mau") for s in row):
        header_row = idx
        break

if header_row is not None:
    df = pd.read_excel(file_path, header=header_row, engine="openpyxl")
else:
    # fallback: assume ml/mAU headers on row 3, data starts on row 4
    df = pd.read_excel(file_path, header=2, engine="openpyxl")

# === IDENTIFY COLUMNS ===
cols_lower = [c.lower().strip() for c in df.columns]
ml_col = None
mau_col = None
for i, cname in enumerate(cols_lower):
    if "ml" in cname:
        ml_col = df.columns[i]
    if "mau" in cname:
        mau_col = df.columns[i]

if ml_col is None or mau_col is None:
    print("⚠️ Could not detect headers 'ml' and 'mAU' automatically, trying fallback.")
    ml_col, mau_col = df.columns[0], df.columns[1]

# === CLEAN AND PREPARE DATA ===
df[ml_col] = pd.to_numeric(df[ml_col], errors="coerce")
df[mau_col] = pd.to_numeric(df[mau_col], errors="coerce")
df_clean = df[[ml_col, mau_col]].dropna()
df_clean.columns = ["ml", "mAU"]
df_clean = df_clean.sort_values("ml").reset_index(drop=True)

print(f"Detected columns: {ml_col} (x), {mau_col} (y)")
print(f"Data rows used: {len(df_clean)}")
print(f"ml range: {df_clean['ml'].min()} — {df_clean['ml'].max()}")
print(f"mAU range: {df_clean['mAU'].min()} — {df_clean['mAU'].max()}")

# === PLOT ===
plt.figure(figsize=(6, 4), facecolor='none')
plt.plot(df_clean["ml"], df_clean["mAU"], linewidth=1.2, color="black")
plt.xlabel("Retention Volume (mL)")
plt.ylabel("mAU")
plt.title("Chromatogram")
plt.tight_layout()
plt.savefig("chromatogram.png", transparent=True, bbox_inches='tight', dpi=300)
plt.show()