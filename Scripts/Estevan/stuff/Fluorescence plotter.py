import pandas as pd
import matplotlib.pyplot as plt

# --- User inputs ---
filename = r"D:\UTAustin Summer 2025\xlsx files\mblueberry2.xlsx"

# --- Read Excel file ---
# skip the first row (A1/B1), and let pandas read headers from row 2
data = pd.read_excel(filename, skiprows=1)

# --- Quick check ---
print(data.head())

# --- Plot ---
plt.figure(figsize=(8, 5))
plt.plot(data.iloc[:, 0], data.iloc[:, 1], linewidth=2)

# --- Labels and formatting ---
plt.xlabel("Wavelength (nm)", fontsize=12)
plt.ylabel("Intensity (A.U.)", fontsize=12)
plt.title("Fluorescence Spectrum", fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.show()
