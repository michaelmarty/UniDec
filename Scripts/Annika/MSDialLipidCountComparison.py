
import pandas as pd

# ------------------------------
# User files
# ------------------------------
uncorrected_file = r"Z:\Group Share\Annika\IQ-X\UntDDA Tests HEKpos\Peak Masters\NegPeakMaster_HEKneg_IT_IT_NoCorr.csv"
corrected_file   = r"Z:\Group Share\Annika\IQ-X\UntDDA Tests HEKpos\New Libs Tests\Peak Masters\NegPeakMaster_HEKneg_ITIT_WithCorr.csv"

# ------------------------------
# Load only Metabolite name
# ------------------------------
df_uncorr = pd.read_csv(uncorrected_file, usecols=["Metabolite name"], engine='python', on_bad_lines='skip')
df_corr   = pd.read_csv(corrected_file, usecols=["Metabolite name"], engine='python', on_bad_lines='skip')

# ------------------------------
# Load color map (optional, used to label "Other")
# ------------------------------
colormap_file = r"Z:\\Group Share\\Annika\\Lipid_ColorMap.xlsx"
color_df = pd.read_excel(colormap_file)
tmp_class = color_df["Class"].astype(str).str.strip()
tmp_color = color_df["Color"].astype(str).str.strip()
color_map = dict(zip(tmp_class, tmp_color))
def_color = "#808080"

# ------------------------------
# Class extraction (use your full plotting function)
# ------------------------------
def get_class(name):
    if not isinstance(name, str):
        return "Other"
    name = name.strip()

    if name.startswith("Cer "):
        return "Cer"
    if name.startswith("SM "):
        return "SM"
    if name.startswith(("AHexCer", "Hex2Cer", "Hex3Cer", "SHexCer")):
        return "HexCer"
    if name.startswith("GM3"):
        return "HexCer"
    if name.startswith(("PI-Cer", "PE-Cer")):
        return "Cer"
    if name.startswith(("DGCC", "DGDG", "DGGA", "ADGGA", "MGDG", "SQDG")):
        return "DG"
    if name.startswith(("DLCL", "MLCL", "HBMP")):
        return "CL"
    if "d7" in name.lower():
        return "Standard"
    if name.startswith("CL ") and "14:0_14:0_14:0_14:0" in name.lower():
        return "Standard"
    if name.startswith("CoQ"):
        return "CoQ"
    if name.startswith("LPS-N"):
        return "LPS"
    if name.startswith(("PIM", "Ac2PIM", "Ac4PIM", "SPB", "SL", "BA")):
        return "Other"
    for sep in ["(", " "]:
        if sep in name:
            return name.split(sep)[0]
    return name

# Apply class extraction
df_uncorr["Class"] = df_uncorr["Metabolite name"].apply(get_class).str.strip()
df_corr["Class"]   = df_corr["Metabolite name"].apply(get_class).str.strip()

# Apply color map fallback (to match plotting logic)
df_uncorr["Class"] = [c if c in color_map else "Other" for c in df_uncorr["Class"]]
df_corr["Class"]   = [c if c in color_map else "Other" for c in df_corr["Class"]]

# ------------------------------
# Count lipids per class
# ------------------------------
count_uncorr = df_uncorr["Class"].value_counts()
count_corr   = df_corr["Class"].value_counts()

all_classes = set(count_uncorr.index).union(count_corr.index)

# ------------------------------
# Compare and print differences
# ------------------------------
print("Comparison of lipid classes: uncorrected vs corrected\n")
for cls in sorted(all_classes):
    n_uncorr = count_uncorr.get(cls, 0)
    n_corr   = count_corr.get(cls, 0)
    diff = n_corr - n_uncorr  # corrected minus uncorrected
    if diff > 0:
        print(f"{diff} more {cls} lipids in corrected than uncorrected")
    elif diff < 0:
        print(f"{-diff} more {cls} lipids in uncorrected than corrected")
    else:
        print(f"Same number of {cls} lipids in both files ({n_uncorr})")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------
# Example: using the counts we already computed
# ------------------------------
# Assuming you already have:
# count_uncorr = df_uncorr["Class"].value_counts()
# count_corr   = df_corr["Class"].value_counts()

all_classes = sorted(set(count_uncorr.index).union(count_corr.index))
uncorr_counts = [count_uncorr.get(cls, 0) for cls in all_classes]
corr_counts   = [count_corr.get(cls, 0) for cls in all_classes]

# ------------------------------
# Print totals
# ------------------------------
total_uncorr = sum(uncorr_counts)
total_corr   = sum(corr_counts)
print(f"Total lipids identified in uncorrected file: {total_uncorr}")
print(f"Total lipids identified in corrected file:   {total_corr}")
# ------------------------------
# Create stacked bar chart: uncorrected vs corrected
# ------------------------------
x = np.arange(len(all_classes))  # positions for each class
width = 0.35  # bar width

fig, ax = plt.subplots(figsize=(12, 6))

# Bars for uncorrected and corrected
bars1 = ax.bar(x - width/2, uncorr_counts, width, label='Uncorrected', color='skyblue')
bars2 = ax.bar(x + width/2, corr_counts, width, label='Corrected', color='orange')

# Labels
ax.set_xticks(x)
ax.set_xticklabels(all_classes, rotation=45, ha='right')
ax.set_ylabel("Number of lipids")
ax.set_title("Lipid counts per class: Uncorrected vs Corrected in IQX HEKneg IT_IT")
ax.legend()

# Add values on top of bars
for bar in bars1 + bars2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 0.1, f'{int(height)}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.show()
