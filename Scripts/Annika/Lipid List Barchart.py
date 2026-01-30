import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# User input
# =========================
csv_path = "Z:\\Group Share\\Annika\\Stellar\\HEK\\tMS2\\Official Library\\HEK_Neg_Full_tMS2.csv"
color_map_path = "Z:\\Group Share\\Annika\\Lipid_ColorMap.xlsx"
threshold = 1.5  # percent cutoff for "Other"

# =========================
# Load data
# =========================
df = pd.read_csv(csv_path)
df.columns = df.columns.str.strip()

colors_df = pd.read_excel(color_map_path)
colors_df.columns = colors_df.columns.str.strip()

# Expect columns: Lipid_Class | Color
color_dict = dict(zip(colors_df["Lipid Class"], colors_df["Color"]))

# =========================
# Extract lipid class
# =========================
df["Lipid Class"] = df["Compound"].str.split().str[0]

# =========================
# Count lipids per class
# =========================
class_counts = df["Lipid Class"].value_counts()
class_percent = (class_counts / class_counts.sum()) * 100

# =========================
# Collapse rare classes
# =========================
major_percent = class_percent[class_percent >= threshold]
other_percent = class_percent[class_percent < threshold].sum()

major_counts = class_counts[class_percent >= threshold]
other_count = class_counts[class_percent < threshold].sum()

# Sort main classes largest → smallest
major_percent = major_percent.sort_values(ascending=False)
major_counts = major_counts.loc[major_percent.index]

# Append "Other"
if other_percent > 0:
    major_percent = pd.concat([major_percent, pd.Series({"Other": other_percent})])
    major_counts = pd.concat([major_counts, pd.Series({"Other": other_count})])

# Reverse for barh layout (largest top, Other bottom)
major_percent = major_percent[::-1]
major_counts = major_counts.loc[major_percent.index]

# =========================
# Assign colors
# =========================
missing = set(major_percent.index) - set(color_dict.keys())
if missing:
    print(f"⚠️ Missing colors for: {missing}")

bar_colors = [
    color_dict.get(lipid, "#BFBFBF")
    for lipid in major_percent.index
]

# =========================
# Plot
# =========================
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14

plt.figure(figsize=(8, 6))
bars = plt.barh(major_percent.index, major_percent.values, color=bar_colors)

plt.xlabel("Percentage of Total Lipids (%)")
plt.ylabel("Lipid Class")
plt.title("Lipid Class Composition (by Count)")

# =========================
# Add absolute counts on right
# =========================
xmax = major_percent.max()

for bar, lipid in zip(bars, major_percent.index):
    count = int(major_counts[lipid])
    x_pos = bar.get_width() + xmax * 0.01  # small offset
    y_pos = bar.get_y() + bar.get_height() / 2
    plt.text(x_pos, y_pos, f"{count}", va="center", ha="left", fontsize=12)

plt.tight_layout()

# plt.savefig("Lipid_Class_Percentage_ByCount.pdf")
plt.show()