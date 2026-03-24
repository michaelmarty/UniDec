import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14
# ---- FILE PATHS ----
chol_files = [
    r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\369_20260313_AS_04_CB2_FTdiscs_tMS2_APCI.csv",
    r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\369_20260313_AS_02_CB2_Ediscs_tMS2_APCI.csv"
]

is_files = [
    r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\376_20260313_AS_04_CB2_FTdiscs_tMS2_APCI.csv",
    r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\376_20260313_AS_02_CB2_Ediscs_tMS2_APCI.csv"
]

# ---- USER SETTINGS ----
chol_color = "cornflowerblue"
is_color = "lightsteelblue"
is_lw = 2
chol_lw = 2.5
is_lw = 2.5
titles = ["CB2 FT Sample", "CB2 E Sample"]


plt.figure(figsize=(5, 8))

# Loop in reverse for top-bottom flip
for i, (chol_f, is_f, title) in enumerate(zip(chol_files[::-1], is_files[::-1], titles[::-1]), start=1):
    plt.subplot(len(chol_files), 1, i)
    plt.subplots_adjust(hspace=2)

    plt.legend()
    # Cholesterol
    df_chol = pd.read_csv(chol_f)
    plt.plot(df_chol["Time(min)"], df_chol["Relative Abundance"],
             color=chol_color, linewidth=chol_lw, label="Cholesterol")

    # Internal Standard
    df_is = pd.read_csv(is_f)
    plt.plot(df_is["Time(min)"], df_is["Relative Abundance"],
             color=is_color, linewidth=is_lw, label="Internal Standard")
    plt.ylim(bottom=0)
    plt.ylabel("Relative Abundance")
    plt.title(title)
    plt.xlim(7,8.5)
    # Scientific notation without using ax
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.legend(loc="upper right")
plt.xlabel("Time (min)")
plt.tight_layout()
plt.show()