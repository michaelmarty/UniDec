import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import numpy as np

# ------------------------------
# User files
# ------------------------------
file = "Z:\Group Share\\Annika\\Stellar\\Untargeted DDA\\HEK\\New Libs Tests\\NegPeakMaster_CorrLib_Ref.csv"
colormap_file = r"Z:\\Group Share\\Annika\\Lipid_ColorMap.xlsx"

# ------------------------------
# Load CSV
# ------------------------------
df = pd.read_csv(file, delimiter=",")

blacklist = "RIKEN|Unknown|low|no MS2"
df = df[~df["Metabolite name"].astype(str).str.contains(blacklist, case=False, na=False)]

# Remove ALL fatty acids (FA, FFA, Fatty acid)
fa_pattern = r"^FA\b|^FFA\b|fatty acid"
df = df[~df["Metabolite name"].astype(str).str.contains(fa_pattern, case=False, na=False)]
# ------------------------------
# Load color map
# ------------------------------
color_df = pd.read_excel(colormap_file)
tmp_class = color_df["Lipid Class"].astype(str).str.strip()
tmp_color = color_df["Color"].astype(str).str.strip()
color_map = dict(zip(tmp_class, tmp_color))

# ------------------------------
# Class extraction
# ------------------------------

def get_class(name):
    if not isinstance(name, str):
        return "Other"
    name = name.strip()

    # -------------------
    # Straightforward sphingolipids
    # -------------------
    if name.startswith("Cer "):
        return "Cer"
    if name.startswith("SM "):
        return "SM"

    # HexCer variants
    if name.startswith(("AHexCer", "Hex2Cer", "Hex3Cer", "SHexCer")):
        return "HexCer"

    # Gangliosides → HexCer
    if name.startswith("GM3"):
        return "HexCer"

    # PE-Cer, PI-Cer → Cer
    if name.startswith(("PI-Cer", "PE-Cer")):
        return "Cer"

    # -------------------
    # Glycerolipids: DG cluster
    # -------------------
    if name.startswith(("DGCC", "DGDG", "DGGA", "ADGGA", "MGDG", "SQDG")):
        return "DG"

    # -------------------
    # Cardiolipin family → CL
    # -------------------
    if name.startswith(("DLCL", "MLCL", "HBMP")):
        return "CL"

    if "d7" in name.lower():
        return "Standard"

    if name.startswith("CL ") and "14:0_14:0_14:0_14:0" in name.lower():
        return "Standard"
    # -------------------
    # CoQ stays as CoQ
    # -------------------
    if name.startswith("CoQ"):
        return "CoQ"
    if name.startswith("LPS-N"):
        return "LPS"


    # -------------------
    # Bacterial-only classes stay in Other
    # -------------------
    if name.startswith(("PIM", "Ac2PIM", "Ac4PIM", "SPB", "SL", "BA")):
        return "Other"

    # -------------------
    # Default assignment (first token before space or parenthesis)
    # -------------------
    for sep in ["(", " "]:
        if sep in name:
            return name.split(sep)[0]

    return name

# Apply extraction
df["Class"] = df["Metabolite name"].apply(get_class).str.strip()

# ------------------------------
# Assign final class/color
# ------------------------------
def_color = "#808080"
final_classes = []
colors = []

for c in df["Class"]:
    if c in color_map:
        final_classes.append(c)
        colors.append(color_map[c])
    else:
        final_classes.append("Other")
        colors.append(def_color)

df["Class"] = final_classes
df["Color"] = colors
# ------------------------------
# Count lipids per class
# ------------------------------
class_counts = df["Class"].value_counts()
print("\nNumber of lipids per class:")
for cls, count in class_counts.items():
    print(f" • {cls}: {count}")
# ------------------------------
# Diagnostics
# ------------------------------
print("Unique extracted classes:", sorted(df["Class"].unique()))
print("Classes in color map:", sorted(color_map.keys()))
print("Number of points plotted:", len(df))

# List "Other"
other_df = df[df["Class"] == "Other"]
print("\nLipids labeled as 'Other' (unique names):")
for name in sorted(other_df["Metabolite name"].unique()):
    print(" •", name)
print(f"\nTotal unique 'Other' lipids: {len(other_df['Metabolite name'].unique())}")

# ------------------------------
# Plot setup
# ------------------------------
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 18

# Legend order: alphabetical, but put Other last if present
unique_classes = sorted([c for c in df["Class"].unique() if c != "Other"])
if "Other" in df["Class"].unique():
    unique_classes.append("Other")

# Build legend
custom_legend = [
    Line2D([0], [0], marker='o', color=color_map.get(cls, def_color), label=cls,
           markerfacecolor=color_map.get(cls, def_color), linestyle='None', markersize=10)
    for cls in unique_classes
]

# ------------------------------
# Plot
# ------------------------------
plt.figure(figsize=(12, 8))
plt.scatter(df["Average Rt(min)"], df["Average Mz"], c=df["Color"], s=80)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Retention Time (min)")
plt.ylabel("m/z")
plt.title("Stellar UntDDA ESI-")
plt.xticks(np.arange(0.0, 18, 2.5))
plt.yticks(np.arange(0.0, 1600.1, 200))

# Add total number of lipids on the plot
total_lipids = len(df)
plt.text(
    0.95, 0.95, f"Total lipids: {total_lipids}",
    horizontalalignment='right',
    verticalalignment='top',
    transform=plt.gca().transAxes,
    fontsize=14,
    bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=5)
)

ncol = min(6, len(unique_classes))
plt.legend(handles=custom_legend, loc='upper center', bbox_to_anchor=(0.5, -0.08),
           ncol=ncol, frameon=False)

plt.subplots_adjust(bottom=0.25)
plt.tight_layout()
plt.show()
