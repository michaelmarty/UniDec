import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# LOAD
# =========================

file = r"Z:\Group Share\Annika\Stellar\TAMU\Natural\NaturalSamples_Molecule Transition Results.csv"
df = pd.read_csv(file)
df.columns = df.columns.str.strip()
control_sample = "BPE"
samples_to_plot = ["1 N", "2 N", "3 N", "4 N", "5 N"]

# REQUIRED: class → internal standard name mapping
internal_standard_map = {
    "Cer": "IS Cer",
    "CL": "IS CL",
    "DG": "IS DG",
    "HexCer_NS": "IS Cer",
    "PC": "IS PC",
    "PE": "IS PE",
    "PG": "IS PG",
    "PI": "IS PI",
    "PS": "IS PS",
    "SM": "IS SM",
    "TG": "IS TG",
    "LPC": "IS LPC",
    "LPE": "IS LPE",
    "LPG": "IS LPG",
    "LPI": "IS LPI",
    "LPS": "IS LPS",
    "PA": "IS PA",
    "Cer_NS": "IS Cer",
    "EtherLPC": "IS LPC",
    "EtherPC": "IS PC",
    "EtherPE": "IS PE",
    "Hex2Cer": "IS Cer",
    "HexCer_HS": "IS Cer",
    "HexCer_NDS": "IS Cer",
    "HexCer_NS": "IS Cer",
    "SHexCer": "IS Cer",
}

# -----------------------------
# LOAD
# -----------------------------


df.columns = df.columns.str.strip()
df["Replicate Name"] = df["Replicate Name"].astype(str).str.strip()
df["Molecule List Name"] = df["Molecule List Name"].astype(str).str.strip()
df["Molecule"] = df["Molecule"].astype(str).str.strip()

# -----------------------------
# STEP 1: collapse fragment ions
# -----------------------------
df = df.groupby(
    ["Molecule", "Molecule List Name", "Replicate Name"],
    as_index=False
)["Area"].sum()

# -----------------------------
# STEP 2: FORCE class → IS assignment (NO GUESSING)
# -----------------------------

df["IS_class"] = df["Molecule List Name"].map(internal_standard_map)

# keep ONLY analytes that have a defined IS mapping
an_df = df[df["IS_class"].notna()].copy()

# IS rows are anything starting with IS
is_df = df[df["Molecule"].astype(str).str.startswith("IS")].copy()

print("Internal standard rows found:", len(is_df))

# -----------------------------
# STEP 3: build IS lookup (FORCED MATCH)
# -----------------------------
is_lookup = []

for cls, is_name in internal_standard_map.items():

    is_subset = is_df[
        is_df["Molecule"].str.contains(is_name.replace("IS ", ""), na=False)
    ]

    tmp = is_subset.groupby(["Replicate Name"])["Area"].mean().reset_index()
    tmp["Molecule List Name"] = cls
    tmp = tmp.rename(columns={"Area": "IS_area"})

    is_lookup.append(tmp)

is_lookup = pd.concat(is_lookup, ignore_index=True)

# -----------------------------
# STEP 4: merge IS onto analytes
# -----------------------------
an_df = an_df.merge(
    is_lookup,
    on=["Molecule List Name", "Replicate Name"],
    how="left"
)
# -----------------------------
# STEP 5: normalize
# -----------------------------
an_df["normalized"] = an_df["Area"] / an_df["IS_area"]

# drop broken rows
an_df = an_df.dropna(subset=["normalized"])

df_mat = an_df.copy()
df_mat["log2_intensity"] = np.log2(df_mat["normalized"] + 1e-12)

heatmap_df = df_mat.pivot_table(
    index="Molecule",
    columns="Replicate Name",
    values="log2_intensity",
    aggfunc="mean"
)

plt.figure(figsize=(10, 8))
plt.imshow(heatmap_df.fillna(0), aspect="auto", cmap="viridis")
plt.colorbar(label="log2(IS-normalized intensity)")

plt.yticks(range(len(heatmap_df.index)), heatmap_df.index, fontsize=6)
plt.xticks(range(len(heatmap_df.columns)), heatmap_df.columns, rotation=45)

plt.title("Lipid Abundance Heatmap (log2 IS-normalized)")
plt.tight_layout()
plt.show()

class_df = df_mat.groupby(
    ["Molecule List Name", "Replicate Name"]
)["normalized"].sum().reset_index()

class_pivot = class_df.pivot(
    index="Molecule List Name",
    columns="Replicate Name",
    values="normalized"
).fillna(0)

class_pivot.plot(kind="bar", figsize=(12,6))

plt.ylabel("Total IS-normalized abundance")
plt.title("Lipid Class Composition Across Samples")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

control = "BPE"

delta = heatmap_df.subtract(heatmap_df[control], axis=0)

plt.figure(figsize=(10, 8))
plt.imshow(delta.fillna(0), aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
plt.colorbar(label="Δ log2 intensity vs BPE")

plt.yticks(range(len(delta.index)), delta.index, fontsize=6)
plt.xticks(range(len(delta.columns)), delta.columns, rotation=45)

plt.title("Differential Lipid Abundance (vs BPE)")
plt.tight_layout()
plt.show()

import numpy as np

control = an_df[an_df["Replicate Name"] == "BPE"]
control_mean = control.groupby("Molecule List Name")["normalized"].mean()

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()

for i, sample in enumerate(samples_to_plot):

    ax = axes[i]

    sub = an_df[an_df["Replicate Name"] == sample]

    if sub.empty:
        ax.set_title(sample)
        ax.axis("off")
        continue

    sample_mean = sub.groupby("Molecule List Name")["normalized"].mean()

    fc = (sample_mean / control_mean).replace([np.inf, -np.inf], np.nan)
    log2fc = np.log2(fc)
    intensity = np.log10(sample_mean)

    ax.scatter(log2fc, intensity, s=10)

    ax.axvline(0, color="black", linewidth=1)
    ax.set_title(sample)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("log10 intensity")

plt.tight_layout()
plt.show()