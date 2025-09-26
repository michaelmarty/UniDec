import os
import pandas as pd
import matplotlib.pyplot as plt

# Set global font to Arial, size 28 for all plot text
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# List of file paths to chromatogram CSV files
file_paths = [
    r'Z:\Group Share\BHT\Akta FPLC Data\20250916\20250916_BHT_HIC_Butyl S-FF_12uM MSP D1T0(-)_500 uL 001.csv',
    r'Z:\Group Share\BHT\Akta FPLC Data\20250917\20250917_BHT_HIC_Buytl S-FF_15uMNanodisc 001.csv',
    r'Z:\Group Share\BHT\Akta FPLC Data\20250917\20250917_BHT_HIC_Butyl S-FF_MPBlank 001.csv'
]

# Labels for each chromatogram trace
labels = [
    '12uM MSP D1T0(-)',
    '15uM Nanodisc',
    'MP Blank'
]

# Colors for each trace
colors = ['blue', 'green', 'orange']

# Read all CSV files, skipping the first 2 rows, using UTF-16-LE encoding and tab separator
dfs = [pd.read_csv(fp, skiprows=2, encoding="utf-16-le", sep='\t') for fp in file_paths]

# Extract the first two columns (Volume, UV) from each dataframe for plotting
UV_datas = [df.iloc[:, :2] for df in dfs]

# Extract columns 4 and 5 (%B data) from the first file only for the secondary y-axis
Concentration_B_Data_1 = dfs[0].iloc[:, 4:6]

# Create a figure and primary axis with 16:9 aspect ratio
fig, ax1 = plt.subplots(figsize=(16, 9))

# Plot UV data for each file on the primary y-axis
for UV_Data, label, color in zip(UV_datas, labels, colors):
    ax1.plot(UV_Data.iloc[:, 0], UV_Data.iloc[:, 1], color=color, label=label)

# Set axis labels and limits for the primary y-axis
ax1.set_xlabel('Volume (mL)')
ax1.set_ylabel('UV Absorbance (mAU)')
ax1.tick_params(axis='y')
ax1.set_xlim(left=0, right=40)
ax1.set_ylim(bottom=-2)
ax1.legend(loc='upper left', fontsize=18)

# Create a secondary y-axis for %B data from the first file
ax2 = ax1.twinx()
ax2.plot(Concentration_B_Data_1.iloc[:, 0], Concentration_B_Data_1.iloc[:, 1], color='red', label='% B')
ax2.set_ylabel('% Mobile Phase B', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(left=0, right=40)

# Save the plot as a PNG in the same directory as the first CSV file
output_dir = os.path.dirname(file_paths[0])
output_path = os.path.join(output_dir, 'Akta_Chromatogram_Comparison.png')
plt.savefig(output_path, bbox_inches='tight', dpi=300)
plt.close()