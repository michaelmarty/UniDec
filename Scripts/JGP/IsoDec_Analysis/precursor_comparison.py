import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Set backend to Agg
mpl.use('WxAgg')

classic_file = "C:\\Users\\johnn\\Desktop\\MartyLab\\20170105_L_MaD_ColC4_Ecoli20161116-12_01_precursors_classic.txt"
isodec_file = "C:\\Users\\johnn\\Desktop\\MartyLab\\20170105_L_MaD_ColC4_Ecoli20161116-12_01_precursors_isodec.txt"

classic_df = pd.read_csv(classic_file, sep='\t')
isodec_df = pd.read_csv(isodec_file, sep='\t')

classic_masses = classic_df["Precursor Mass (Da)"]
isodec_masses = isodec_df["Precursor Mass (Da)"]

classic_charges = classic_df["Precursor Charge"]
isodec_charges = isodec_df["Precursor Charge"]

# classic_charges_over3_bool = classic_charges > 3
# isodec_charges_over3_bool = isodec_charges > 3
#
# classic_masses = classic_masses[classic_charges_over3_bool]
# isodec_masses = isodec_masses[isodec_charges_over3_bool]

classic_over20k = 0
isodec_over20k = 0

#manually make a histogram of the proteoform masses
min_mass = 0
max_mass = 30000
bin_width = 1000
bins = np.arange(min_mass, max_mass, bin_width)
file1_hist = np.zeros(len(bins) - 1)
file2_hist = np.zeros(len(bins) - 1)

for mass in classic_masses:
    if mass > 20000:
        classic_over20k += 1
    for j in range(len(bins) - 1):
        if mass >= bins[j] and mass < bins[j + 1]:
            file1_hist[j] += 1
            break

for mass in isodec_masses:
    if mass > 20000:
        isodec_over20k += 1
    for j in range(len(bins) - 1):
        if mass >= bins[j] and mass < bins[j + 1]:
            file2_hist[j] += 1
            break

print("Classic > 20 kDa: ", classic_over20k)
print("IsoDec > 20 kDa: ", isodec_over20k)

#Normalize counts to the total number of PSMs
# file1_hist /= len(classic_masses)
# file2_hist /= len(isodec_masses)

bar_width = 0.35

index = np.arange(len(bins[:-1]))

fig, ax = plt.subplots()

bar1 = ax.bar(index, file1_hist, bar_width, label="Classic")
bar2 = ax.bar(index + bar_width, file2_hist, bar_width, label="IsoDec")

ax.set_xlabel("Proteoform Mass (Da)")
ax.set_ylabel("Count")

ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels((str(bins[i]) for i in range(len(bins[:-1]))), rotation=45, ha='right')

ax.legend(fontsize=32)
plt.show()
exit()

#manually make a histogram of the observed charges
charges = np.arange(1,50)
classic_counts = np.zeros(len(charges))
isodec_counts = np.zeros(len(charges))

for charge in classic_charges:
    for i in range(len(charges)):
        if charge == charges[i]:
            classic_counts[i] += 1

for charge in isodec_charges:
    for i in range(len(charges)):
        if charge == charges[i]:
            isodec_counts[i] += 1

bar_width = 0.35

index = np.arange(len(charges))

fig, ax = plt.subplots()

bar1 = ax.bar(index, classic_counts, bar_width, label="Classic")
bar2 = ax.bar(index + bar_width, isodec_counts, bar_width, label="IsoDec")

ax.set_xlabel("Charge")
ax.set_ylabel("Count")

ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(charges)

ax.legend()
plt.show()

