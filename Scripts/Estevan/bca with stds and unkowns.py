import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# -------------------------------------------------------------
# 🧠 USER INPUT SECTION
file_path = r"D:\bcaallmsp.xlsx"   # ← your Excel file path
sheet_name = 0  # change if needed
# -------------------------------------------------------------

# Load Excel (no headers)
df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)

# -------------------------------------------------------------
# 🧪 Standards (positions based on Excel sheet)
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

# -------------------------------------------------------------
# 🧩 Unknowns (positions + dilutions)
unknowns = {
    "U1": {
        "1:10":   ['K2','K3','K4'],
        "1:100":  ['L2','L3','L4'],
        "1:1000": ['M2','M3','M4']
    },
    "U2": {
        "1:10":   ['B5','B6','B7'],
        "1:100":  ['C5','C6','C7'],
        "1:1000": ['D5','D6','D7']
    },
    "U3": {
        "1:10":   ['E5','E6','E7'],
        "1:100":  ['F5','F6','F7'],
        "1:1000": ['G5','G6','G7']
    },
    "U4": {
        "1:10":   ['H5','H6','H7'],
        "1:100":  ['I5','I6','I7'],
        "1:1000": ['J5','J6','J7']
    },
    "U5": {
        "1:10":   ['K5','K6','K7'],
        "1:100":  ['L5','L6','L7'],
        "1:1000": ['M5','M6','M7']
    }
}

# -------------------------------------------------------------
# 🧮 Helper: Excel coordinate → DataFrame index
def excel_to_idx(cell):
    col_letters = ''.join([c for c in cell if c.isalpha()])
    row_numbers = ''.join([c for c in cell if c.isdigit()])
    col_idx = sum([(ord(c.upper()) - 64) * (26 ** i) for i, c in enumerate(col_letters[::-1])])
    row_idx = int(row_numbers)
    return row_idx - 1, col_idx - 1  # zero-based for pandas

# -------------------------------------------------------------
# 🧪 Extract standards
conc, means, sds = [], [], []
for c, cells in standards.items():
    vals = []
    for cell in cells:
        r, col = excel_to_idx(cell)
        vals.append(float(df.iat[r, col]))
    vals = np.array(vals)
    conc.append(c)
    means.append(np.mean(vals))
    sds.append(np.std(vals, ddof=1))

conc = np.array(conc)
means = np.array(means)
sds = np.array(sds)

# -------------------------------------------------------------
# 🧪 Blank correction
blank_mean = means[conc == 0][0]
means_corrected = means - blank_mean
print(f"\nBlank (0 µg/mL) mean absorbance = {blank_mean:.4f}")

# -------------------------------------------------------------
# 📈 Linear regression (omit blank)
mask = conc > 0
slope, intercept, r_value, p_value, std_err = linregress(conc[mask], means_corrected[mask])
r2 = r_value ** 2

# Fitted line for plotting
x_fit = np.linspace(0, max(conc)*1.1, 200)
y_fit = intercept + slope * x_fit

# -------------------------------------------------------------
# 🎨 Plot styling
plt.rcParams.update({'font.family': 'Arial', 'font.size': 10})

plt.figure(figsize=(7,5))
plt.errorbar(conc, means_corrected, yerr=sds, fmt='o', capsize=4, color='black', label='Standards (mean ± SD)')
plt.plot(x_fit, y_fit, '--', color='black', linewidth=1, label=f'Fit: y = {slope:.4f}x + {intercept:.4f}\n$R^2$ = {r2:.4f}')

# -------------------------------------------------------------
# 🧪 Process unknowns
unk_results = []
colors = plt.cm.tab10.colors
for i, (unk_name, dil_dict) in enumerate(unknowns.items()):
    color = colors[i % len(colors)]
    for dil, cells in dil_dict.items():
        vals = []
        for cell in cells:
            r, col = excel_to_idx(cell)
            vals.append(float(df.iat[r, col]))
        vals = np.array(vals)
        mean_abs = np.mean(vals) - blank_mean  # blank corrected
        sd_abs = np.std(vals, ddof=1)
        conc_measured = (mean_abs - intercept) / slope  # µg/mL in the *diluted* sample

        # Get dilution factor as integer
        dil_factor = int(dil.split(':')[1])
        conc_original = conc_measured * dil_factor  # back-calculate original sample conc

        unk_results.append([unk_name, dil, mean_abs, sd_abs, conc_measured, conc_original])

        # Plot
        plt.errorbar(conc_measured, mean_abs, yerr=sd_abs, fmt='s', markersize=6,
                     color=color, label=f'{unk_name} {dil}')

# -------------------------------------------------------------
# 🎨 Final plot formatting
plt.xlabel('Protein concentration (µg/mL)')
plt.ylabel('Absorbance (562 nm, blank-corrected)')
plt.title('BCA Standard Curve with Unknowns', pad=10)
plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# -------------------------------------------------------------
# 📊 Output table
unk_df = pd.DataFrame(
    unk_results,
    columns=["Unknown", "Dilution", "Mean Abs (corr)", "SD", "Measured Conc (µg/mL)", "Back-calculated Conc (µg/mL)"]
)

print("\n--- Calculated Unknown Concentrations ---")
print(unk_df.to_string(index=False))

# Optional: save to CSV
output_path = file_path.replace('.xlsx', '_BCA_results.csv')
unk_df.to_csv(output_path, index=False)
print(f"\nResults saved to: {output_path}")

# -------------------------------------------------------------
# Regression info
print(f"\nSlope: {slope:.6f}")
print(f"Intercept: {intercept:.6f}")
print(f"R²: {r2:.4f}")
