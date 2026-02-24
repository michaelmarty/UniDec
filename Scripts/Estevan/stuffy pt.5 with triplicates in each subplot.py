import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# === FONT SETTINGS ===
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 10

# === INPUT FILES (triplicates per condition) ===
file_groups = [
    [
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1T0minus_1to80_DMPC 001.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 003D1T0minus_1to80_DMPC 001.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\D1T0(minus)1to80DMPC.xlsx"
    ],
    [
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 003D1T0minus_1to80_DMPC 001.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1t65_DMPC.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\D1sfGFP_422uM_1to60_DMPC.xlsx"
    ],
    [
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 003sfCherryshort_1to80_DMPC_rep2 002.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\Sample3_rep2.xlsx",
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_maybeD1sfCherryshort_minustag_1to150DMPC 001.xlsx"
    ],
    [
        r"D:\UTAustin Summer 2025\xlsx files\Sample4_rep1.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample4_rep2.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample4_rep3.xlsx"
    ],
    [
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_D1halo1toDMPC 001.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample5_rep2.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample5_rep3.xlsx"
    ],
    [
        r"D:\UTAustin Summer 2025\xlsx files\NDtry1 002_E3sfGFPminus_1to350_DMPC 001.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample6_rep2.xlsx",
        #r"D:\UTAustin Summer 2025\xlsx files\Sample6_rep3.xlsx"
    ]
]

# === LABELS FOR EACH CONDITION ===
labels = ["Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6"]

# === CUSTOM COLORS ===
colors = [
    "#440154FF",  # purple
    "#414487FF",  # blue
    "#2A788EFF",  # teal
    "#22A884FF",  # green
    "#7AD151FF",  # lime
    "#FDE725FF"  # yellow-green
]


# === FUNCTION TO READ FILES STARTING AT A4 ===
def read_sec_file(file_path):
    df = pd.read_excel(file_path, skiprows=3, engine="openpyxl")  # starts from row 4 (A4/B4)
    df.columns = [str(c).strip().lower() for c in df.columns]
    if "ml" not in df.columns or "mau" not in df.columns:
        # Assume first two columns are ml/mAU if not labeled
        df = df.iloc[:, :2]
        df.columns = ["ml", "mau"]
    df = df.dropna()
    df["mau"] = df["mau"] / df["mau"].max()  # normalize
    return df


# === PLOTTING ===
fig, axes = plt.subplots(len(file_groups), 1, sharex=True, figsize=(6, 8), dpi=300)

for i, (group, label) in enumerate(zip(file_groups, labels)):
    ax = axes[i]

    for j, file_path in enumerate(group):
        if not os.path.exists(file_path):
            print(f"⚠️ File not found: {file_path}")
            continue
        df = read_sec_file(file_path)
        ax.plot(
            df["ml"],
            df["mau"],
            color=colors[i % len(colors)],
            linewidth=1.4,
            alpha=0.7 + 0.1 * j,  # small transparency difference
            label=f"{label} rep{j + 1}"
        )

    # Styling per subplot
    ax.set_xlim(10, 20)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

# === SHARED LABELS ===
axes[-1].set_xlabel("Elution Volume (mL)")
fig.text(0.03, 0.5, "Normalized Absorbance (mAU)", va="center", rotation="vertical")

# === LAYOUT ===
plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.8)
plt.subplots_adjust(hspace=0.37)

# === SAVE & SHOW ===
plt.savefig("SEC_triplot_normalized.png", dpi=300, bbox_inches="tight", transparent=True)
plt.show()

print("✅ SEC triplicate plots saved as 'SEC_triplot_normalized.png'")
