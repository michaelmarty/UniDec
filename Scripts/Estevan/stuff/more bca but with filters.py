import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.lines as mlines
import os

# =========================================
# 🔧 USER INPUT
# =========================================
file_path = r"D:\bcaallmsp.xlsx"  # <<<--- put your file path here

# Molecular masses for each unknown in Daltons
protein_masses = {
    "Unknown 1": 56799,
    "Unknown 2": 48861,
    "Unknown 3": 47746,
    "Unknown 4": 48104,
    "Unknown 5": 48702
}

# Define standards
standards = {
    2000: ['B2','B3','B4'],
    1500: ['C2','C3','C4'],
    1000: ['D2','D3','D4'],
    750:  ['E2','E3','E4'],
    500:  ['F2','F3','F4'],
    250:  ['G2','G3','G4'],
    125:  ['H2','H3','H4'],
    25:   ['I2','I3','I4'],
    0:    ['J2','J3','J4']
}

# Unknowns layout
unknowns = {
    "Unknown 1": {"1:10": ["K2","K3","K4"], "1:100": ["L2","L3","L4"], "1:1000": ["M2","M3","M4"]},
    "Unknown 2": {"1:10": ["B5","B6","B7"], "1:100": ["C5","C6","C7"], "1:1000": ["D5","D6","D7"]},
    "Unknown 3": {"1:10": ["E5","E6","E7"], "1:100": ["F5","F6","F7"], "1:1000": ["G5","G6","G7"]},
    "Unknown 4": {"1:10": ["H5","H6","H7"], "1:100": ["I5","I6","I7"], "1:1000": ["J5","J6","J7"]},
    "Unknown 5": {"1:10": ["K5","K6","K7"], "1:100": ["L5","L6","L7"], "1:1000": ["M5","M6","M7"]}
}

# =========================================
# Helper
def excel_to_idx(cell):
    col = cell[0].upper()
    row = int(cell[1:]) - 1
    col_idx = ord(col) - ord('A')
    return row, col_idx

# Load data
df = pd.read_excel(file_path, header=None)

# =========================================
# Process standards
std_conc, std_abs = [], []
for conc, wells in standards.items():
    vals = [float(df.iat[excel_to_idx(w)[0], excel_to_idx(w)[1]]) for w in wells]
    std_conc.append(conc)
    std_abs.append(np.mean(vals))

std_conc = np.array(std_conc)
std_abs = np.array(std_abs)
blank_mean = np.mean(std_abs[std_conc == 0])
abs_corrected = std_abs - blank_mean

# Select linear portion automatically
mask = std_conc > 0
slope, intercept, r_value, _, _ = linregress(std_conc[mask], abs_corrected[mask])
pred = slope * std_conc[mask] + intercept
residuals = abs(abs_corrected[mask] - pred)
threshold = 0.10 * np.max(abs_corrected)
linear_mask = residuals < threshold

# Refit using only linear points
slope, intercept, r_value, _, _ = linregress(std_conc[mask][linear_mask], abs_corrected[mask][linear_mask])
pred_linear = slope*std_conc[mask][linear_mask] + intercept

# LOD, LOQ, LOL using linear portion
sigma_blank = np.std([float(df.iat[*excel_to_idx(w)]) for w in standards[0]], ddof=1)
LOD = 3*sigma_blank / slope
LOQ = 10*sigma_blank / slope
residuals_linear = abs(abs_corrected[mask][linear_mask] - pred_linear)
resid_threshold = 0.1*np.max(abs_corrected[mask][linear_mask])
LOL = std_conc[mask][linear_mask][np.where(residuals_linear < resid_threshold)[0][-1]]

# =========================================
# Process unknowns
unknown_results = []
valid_list = []

for i, (unk_name, dil_dict) in enumerate(unknowns.items()):
    for dil, cells in dil_dict.items():
        vals = [float(df.iat[excel_to_idx(w)[0], excel_to_idx(w)[1]]) for w in cells]
        mean_abs = np.mean(vals) - blank_mean
        sd_abs = np.std(vals, ddof=1)
        conc_measured = max((mean_abs - intercept)/slope, 0)
        dil_factor = int(dil.split(':')[1])
        conc_original = conc_measured * dil_factor

        # Correct µM calculation (Daltons)
        conc_uM = conc_original * 1e3 / protein_masses[unk_name]

        in_range = (conc_measured >= LOQ) and (conc_measured <= LOL)
        status = "VALID" if in_range else "OUT OF RANGE"
        if in_range:
            valid_list.append(f"{unk_name} ({dil})")

        unknown_results.append({
            "Unknown": unk_name,
            "Dilution": dil,
            "Mean Abs": round(mean_abs,4),
            "SD Abs": round(sd_abs,4),
            "Measured Conc (µg/mL)": round(conc_measured,2),
            "Back-calculated Conc (µg/mL)": round(conc_original,2),
            "Conc (µM)": round(conc_uM,3),
            "Status": status
        })

# =========================================
# Print and save results BEFORE plotting
unk_df = pd.DataFrame(unknown_results)
print("\n=== Unknown Concentrations ===")
print(unk_df.to_string(index=False))

output_path = os.path.join(os.path.dirname(file_path), "BCA_unknown_concentrations.csv")
unk_df.to_csv(output_path, index=False)
print(f"\nSaved table to: {output_path}")

print("\n=== Valid Unknowns (within LOQ–LOL range) ===")
if valid_list:
    for v in valid_list:
        print(" -", v)
else:
    print("No valid unknowns found within linear range.")

# =========================================
# Plot standards and fit
plt.figure(figsize=(8,5))
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.errorbar(std_conc, abs_corrected, fmt='o', color='blue', label='Standards')
xfit = np.linspace(0, np.max(std_conc)*1.1, 200)
plt.plot(xfit, slope*xfit + intercept, '--', color='black', label=f'Linear fit (R²={r_value**2:.4f})')
plt.axvline(LOD, color='gray', linestyle=':', lw=1, label=f'LOD={LOD:.1f}')
plt.axvline(LOQ, color='gray', linestyle='--', lw=1, label=f'LOQ={LOQ:.1f}')
plt.axvline(LOL, color='gray', linestyle='-.', lw=1, label=f'LOL={LOL:.0f}')

# Plot unknowns
colors = plt.cm.tab10.colors
for i, (unk_name, dil_dict) in enumerate(unknowns.items()):
    base_color = np.array(colors[i % len(colors)])
    valid_color = base_color
    invalid_color = base_color * 0.5 + 0.5

    for dil, cells in dil_dict.items():
        sub = unk_df[(unk_df["Unknown"] == unk_name) & (unk_df["Dilution"] == dil)]
        mean_abs, sd_abs = sub["Mean Abs"].values[0], sub["SD Abs"].values[0]
        conc_measured = sub["Measured Conc (µg/mL)"].values[0]
        status = sub["Status"].values[0]

        color = valid_color if status == "VALID" else invalid_color
        marker = 's' if status == "VALID" else 'o'

        plt.errorbar(conc_measured, mean_abs, yerr=sd_abs, fmt=marker,
                     color=color, markersize=6, capsize=3, label=f"{unk_name} {dil}" if status == "VALID" else "",
                     mec='black' if status == "VALID" else 'none', mew=0.8)

plt.xlabel("Protein Concentration (µg/mL)")
plt.ylabel("Absorbance (blank-corrected)")
plt.title("BCA Standard Curve with Unknowns")
plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show(block=False)  # 👈 non-blocking show
