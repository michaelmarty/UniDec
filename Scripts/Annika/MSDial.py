import pandas as pd
import matplotlib.pyplot as plt

# === User settings ===
input_file = "Z:\\Group Share\\Annika\\Stellar\\EC\\Raw Data\\Book1.xlsx"   # MSDIAL alignment peak table
rt_column = "Average RT (min)"
mz_column = "Average Mz"
intensity_column = "Average peak height"       # or a specific sample column
name_column = "Metabolite name"
try:
    df = pd.read_csv(input_file, sep="\t", encoding="utf-8", engine="python")
except Exception as e:
    print("UTF-8 with tab failed, retrying with latin1...")
    df = pd.read_csv(input_file, sep="\t", encoding="latin1", engine="python")

print("✅ Loaded successfully!")
print("Columns detected:")
print(df.columns.tolist()[:20])  # show first 20 headers
print(df.head())


# === Load data ===
try:
    df = pd.read_csv(input_file, encoding="cp932")
except UnicodeDecodeError:
    df = pd.read_csv(input_file, encoding="latin1")

# Drop rows missing RT or m/z
df = df.dropna(subset=[rt_column, mz_column])

# Remove "Unknown" or "RIKEN" names
if name_column in df.columns:
    df = df[~df[name_column].str.contains("Unknown|RIKEN", case=False, na=False)]

# Normalize intensity for dot scaling
if intensity_column in df.columns:
    df["norm_intensity"] = df[intensity_column] / df[intensity_column].max()
else:
    df["norm_intensity"] = 0.5  # default size if no intensity column

# === Plot ===
plt.figure(figsize=(12, 7))

scatter = plt.scatter(
    df[rt_column],
    df[mz_column],
    s=df["norm_intensity"] * 100,   # dot size
    c=df["norm_intensity"],         # dot color
    cmap="viridis",
    alpha=0.6,
    edgecolors="none"
)

plt.xlabel("Retention Time (min)")
plt.ylabel("m/z")
plt.title("MSDIAL Alignment Dot Plot")

# Add colorbar if intensity was used
if intensity_column in df.columns:
    plt.colorbar(scatter, label="Normalized Intensity")

# Label a subset of metabolites (optional)
if name_column in df.columns:
    for i, row in df.sample(min(20, len(df))).iterrows():  # label up to 20 random metabolites
        plt.text(row[rt_column], row[mz_column], str(row[name_column]),
                 fontsize=7, alpha=0.7)

plt.tight_layout()
plt.show()