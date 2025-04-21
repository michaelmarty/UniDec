import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import matplotlib

matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['pdf.fonttype'] = 42


pdat_path = "C:\\Data\\UPPDemo\\BsAb\\published.txt"
pdat = np.loadtxt(pdat_path)

os.chdir("C:\\Data\\UPPDemo\\BsAb")

results_path = "C:\\Data\\UPPDemo\\BsAb\\BsAb test - Copy_results.xlsx"
results = pd.read_excel(results_path)

rdat = results["BsAb Pairing Calculated (%)"].to_numpy()

# Calculate RMSD
rmsd = np.sqrt(np.mean((rdat - pdat)**2))
print("RMSD: ", rmsd)

residuals = (rdat - pdat) #/ pdat * 100

# set any nan in residuals to 10
residuals[np.isnan(residuals)] = 10

print(np.average(np.abs(residuals)), np.std(np.abs(residuals)))
print(np.average(residuals), np.std(residuals))
print(np.amax(np.abs(residuals)))
b1 = np.abs(residuals) > 3

results["Residuals (%)"] = residuals
results["Published BsAb Pairing (%)"] = pdat

# Write out the samples that are more than 5% different
results[b1].to_excel("outliers.xlsx", index=False)

print(results["Sample name"][b1].values)

# Plot the data
fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

axes[0].plot([20, 100], [20, 100], 'k--')
axes[0].plot(pdat, rdat, 'o')

axes[0].set_xlabel("Published BsAb Pairing (%)")
axes[0].set_ylabel("UPP BsAb Pairing (%)")

# Remove the top and right spines
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
# set lim to 20-100
axes[0].set_xlim(20, 100)
axes[0].set_ylim(20, 100)

# Plot the residuals

axes[1].plot([20, 100], [0, 0], 'k--')
axes[1].plot(pdat, residuals, 'o')
#axes[1].set_xlabel("Published BsAb Pairing (%)")
axes[1].set_ylabel("Residuals")

# Remove the top and right spines
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
# set lim to 20-100
axes[1].set_xlim(20, 100)

plt.tight_layout()
plt.savefig("bsab_compare.pdf")
plt.show()