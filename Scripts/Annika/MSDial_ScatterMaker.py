import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl

# Load your CSV
file = r"Z:\\Group Share\\Annika\\Stellar\\HEK\\20251023_AS_HEK_PosDDA_MSDial.csv"
df = pd.read_csv(file, delimiter=",")


# --- Filter out unwanted names ---
df = df[~df["Name"].astype(str).str.contains("RIKEN|Unknown|low", case=False, na=False)]

# --- Heuristic to extract lipid class from Name ---
def get_class_from_name(name):
    """
    Extracts the core lipid class from a lipid name.
    Groups variants like CoQ8, CoQ9, CoQ10 into 'CoQ'.
    Assumes class is the first part before a parenthesis or space.
    """
    if not isinstance(name, str):
        return "Unknown"

    # Group CoQ variants together
    if name.startswith("CoQ"):
        return "CoQ"

    # You can add more grouping rules here, e.g.:
    # if name.startswith("Cer"): return "Cer"
    # if name.startswith("PC"): return "PC"

    # Default: take text before '(' or space
    for sep in ['(', ' ']:
        if sep in name:
            return name.split(sep)[0]

    return name  # fallback to full name

# --- Apply class extraction ---
df["Class"] = df["Name"].apply(get_class_from_name)

# --- Assign a unique color to each class dynamically ---
unique_classes = df["Class"].unique()
num_classes = len(unique_classes)
cmap = plt.get_cmap("tab20b")  # choose any colormap
class_colors = {cls: cmap(i / max(num_classes-1, 1)) for i, cls in enumerate(unique_classes)}

# --- Prepare plotting data ---
classes = df["Class"].tolist()
rts = df["RT (min)"].tolist()
mzs = df["Precursor m/z"].tolist()

# --- Create legend ---
custom_legend = [
    Line2D([0], [0], marker='o', color=class_colors[cls], label=cls,
           markerfacecolor=class_colors[cls], linestyle='None', markersize=10)
    for cls in unique_classes
]

# --- Plot ---
plt.figure(figsize=(14,10))
for i in range(len(classes)):
    plt.scatter(rts[i], mzs[i], color=class_colors[classes[i]], s=60)  # adjust dot size

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Retention Time (min)", fontsize=14)
plt.ylabel("m/z", fontsize=14)
plt.title("Stellar HEK +ESI", fontsize=16)

# automatic ncol: roughly 6 classes per row
ncol = min(6, len(unique_classes))
plt.legend(handles=custom_legend, loc='upper center', bbox_to_anchor=(0.5, -0.05),
           ncol=ncol, fontsize=20, frameon=False)

# Make space at bottom for the legend
plt.subplots_adjust(bottom=0.35)
plt.tight_layout()
plt.show()
