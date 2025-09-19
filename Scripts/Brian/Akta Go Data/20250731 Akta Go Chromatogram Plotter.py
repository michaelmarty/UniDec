import pandas as pd
import matplotlib.pyplot as plt

# Set global font to Arial, size 28 for all plot text
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Path to the input CSV file (tab-delimited, UTF-16 LE encoding)
file_path = r'Z:\Group Share\BHT\Akta FPLC Data\20250731\20250731_BHT_HIC_2M_ammoniumacetate_ButylFF 001.csv'

# Read the CSV file, skipping the first two rows, using tab as separator
df = pd.read_csv(file_path, skiprows=2, encoding="utf-16-le", sep='\t')

# Extract first two columns for UV data (Volume and mAU)
UV_Data = df.iloc[:, :2]

# Extract fifth and sixth columns for Concentration B data (Volume and Concentration B)
Concentration_B_Data = df.iloc[:, 4:6]

# Create a figure with 16:9 aspect ratio
fig, ax1 = plt.subplots(figsize=(16, 9))

# Plot UV data: Volume (x) vs mAU (y) on primary y-axis
ax1.plot(UV_Data.iloc[:, 0], UV_Data.iloc[:, 1], color='blue', label='UV (mAU)')
ax1.set_xlabel('Volume (mL)')
ax1.set_ylabel('UV Absorbance (mAU)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_xlim(left=0)  # Start x-axis at 0
ax1.set_ylim(bottom=-5)  # Start y-axis at -5 mAU

# Create secondary y-axis for Concentration B
ax2 = ax1.twinx()
# Plot Concentration B: its own Volume (x) vs Concentration B (y)
ax2.plot(Concentration_B_Data.iloc[:, 0], Concentration_B_Data.iloc[:, 1], color='red', label='% B')
ax2.set_ylabel('Concentration B', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(left=0)  # Start x-axis at 0

# Show the plot (no title as requested)
plt.show()