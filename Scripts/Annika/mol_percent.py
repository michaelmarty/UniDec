import pandas as pd
import matplotlib.pyplot as plt


# ============================
# INPUT FILE
# ============================

input_file = r"Z:\Group Share\Annika\Stellar\Cholesterol\CB2 Pellet\Cholesterol_mol%.csv"

df = pd.read_csv(input_file)


# ============================
# SETTINGS
# ============================

starting_sample = "HEK"

cholesterol_name = "Cholesterol"

# edit these to match your lipid names exactly
popc_names = [
    "POPC"
]


# ============================
# CALCULATE LIPID AMOUNT
# ============================

# LC-MS area normalized to internal standard
df["Lipid_umol"] = (
    df["Peak Area"] /
    df["IS Area"]
) * df["IS umol"]


# ============================
# CHOLESTEROL AMOUNT
# ============================

chol = df[
    df["Lipid"].str.contains(
        cholesterol_name,
        case=False,
        na=False
    )
]


chol_summary = (
    chol.groupby("Sample")["Lipid_umol"]
    .sum()
    .reset_index()
)


chol_summary.columns = [
    "Sample",
    "Cholesterol_umol"
]


# ============================
# CHOLESTEROL RECOVERY
# ============================

hek_chol = chol_summary.loc[
    chol_summary["Sample"] == starting_sample,
    "Cholesterol_umol"
].values[0]


chol_summary["Cholesterol_recovery_%"] = (
    chol_summary["Cholesterol_umol"] /
    hek_chol
) * 100


print("\n=== Cholesterol recovery ===")
print(chol_summary)


# ============================
# POPC AMOUNT
# ============================

popc = df[
    df["Lipid"].isin(popc_names)
]


popc_summary = (
    popc.groupby("Sample")["Lipid_umol"]
    .sum()
    .reset_index()
)


popc_summary.columns = [
    "Sample",
    "POPC_umol"
]


# ============================
# CHOLESTEROL : POPC RATIO
# ============================

ratio = pd.merge(
    chol_summary,
    popc_summary,
    on="Sample"
)


ratio["Cholesterol_POPC_ratio"] = (
    ratio["Cholesterol_umol"] /
    ratio["POPC_umol"]
)


print("\n=== Cholesterol : POPC ratio ===")

print(
    ratio[
        [
            "Sample",
            "Cholesterol_POPC_ratio"
        ]
    ]
)


# ============================
# RELATIVE CHOLESTEROL mol%
# ============================

total_lipid = (
    df.groupby("Sample")["Lipid_umol"]
    .sum()
    .reset_index()
)


total_lipid.columns = [
    "Sample",
    "Total_lipid_umol"
]


mol_percent = pd.merge(
    chol_summary,
    total_lipid,
    on="Sample"
)


mol_percent["Relative_Cholesterol_mol%"] = (
    mol_percent["Cholesterol_umol"] /
    mol_percent["Total_lipid_umol"]
) * 100


print("\n=== Relative cholesterol mol% ===")

print(
    mol_percent[
        [
            "Sample",
            "Relative_Cholesterol_mol%"
        ]
    ]
)


# ============================
# PLOTS
# ============================

plt.figure(figsize=(6,4))

plt.bar(
    chol_summary["Sample"],
    chol_summary["Cholesterol_recovery_%"]
)

plt.ylabel("Cholesterol recovery (%)")
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()



plt.figure(figsize=(6,4))

plt.bar(
    ratio["Sample"],
    ratio["Cholesterol_POPC_ratio"]
)

plt.ylabel("Cholesterol / POPC molar ratio")
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()

plt.figure(figsize=(6,4))

plt.bar(
    mol_percent["Sample"],
    mol_percent["Relative_Cholesterol_mol%"]
)

plt.ylabel("Relative cholesterol mol%")
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()