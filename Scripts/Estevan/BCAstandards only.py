import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# -------------------------------------------------------------
# 🧠 USER INPUT SECTION
file_path = r"D:\bcaallmsp.xlsx"  # ← your Excel file path
sheet_name = 0  # change if needed
# -------------------------------------------------------------

# Load Excel (raw values, no headers)
df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)

# Define standards directly by Excel cell positions
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

# Helper: convert Excel coordinates (e.g. 'B2') → numeric indices
def excel_to_idx(cell):
    col_letters = ''.join([c for c in cell if c.isalpha()])
    row_numbers = ''.join([c for c in cell if c.isdigit()])
    col_idx = sum([(ord(c.upper()) - 64) * (26 ** i) for i, c in enumerate(col_letters[::-1])])
    row_idx = int(row_numbers)
    return row_idx - 1, col_idx - 1  # zero-based for pandas

# Extract absorbances
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

# Linear regression
slope, intercept, r_value, p_value, std_err = linregress(conc, means)
r2 = r_value ** 2

# Fit line
x_fit = np.linspace(0, max(conc)*1.1, 200)
y_fit = intercept + slope * x_fit

# -------------------------------------------------------------
# 🎨 Plot styling
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10
})
# -------------------------------------------------------------

plt.figure(figsize=(7,5))
plt.errorbar(conc, means, yerr=sds, fmt='o', capsize=4, color='black', label='Standards (mean ± SD)')
plt.plot(x_fit, y_fit, '--', color='black', linewidth=1, label=f'Fit: y = {slope:.4f}x + {intercept:.4f}\n$R^2$ = {r2:.4f}')
plt.xlabel('Protein concentration (µg/mL)')
plt.ylabel('Absorbance (562 nm)')
plt.title('BCA Standard Curve', pad=10)
plt.legend(frameon=False)
plt.tight_layout()
plt.show()

# Print regression info
print(f"Slope: {slope:.4f}")
print(f"Intercept: {intercept:.4f}")
print(f"R²: {r2:.4f}")
