import os
import pandas as pd
import matplotlib.pyplot as plt

# Set global font to Arial, size 28 for all plot text
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# File paths
file_paths = [
    r'Z:\Group Share\BHT\Akta FPLC Data\20250801\20250801_BHT_HIC_2M_ammoniumacetate_proteinmix_ButylSFF 001.csv',
    r'Z:\Group Share\BHT\Akta FPLC Data\20250801\20250801_BHT_HIC_2M_ammoniumacetate_proteinmix_OctylFF 001.csv'
]
labels = ['UV (ButylSFF)', 'UV (OctylFF)']
colors = ['blue', 'green']

# Read all CSV files
dfs = [pd.read_csv(fp, skiprows=2, encoding="utf-16-le", sep='\t') for fp in file_paths]

# Extract UV data (first two columns) for each file
UV_datas = [df.iloc[:, :2] for df in dfs]
Concentration_B_Data_1 = dfs[0].iloc[:, 4:6]  # Only plot %B from first file

# Create a figure with 16:9 aspect ratio
fig, ax1 = plt.subplots(figsize=(16, 9))

# Plot UV data from all files
for UV_Data, label, color in zip(UV_datas, labels, colors):
    ax1.plot(UV_Data.iloc[:, 0], UV_Data.iloc[:, 1], color=color, label=label)
ax1.set_xlabel('Volume (mL)')
ax1.set_ylabel('UV Absorbance (mAU)')
ax1.tick_params(axis='y')
ax1.set_xlim(left=0, right=50)
ax1.set_ylim(bottom=-5)
ax1.legend(loc='upper right', fontsize=18)

# Create secondary y-axis for Concentration B from first file
ax2 = ax1.twinx()
ax2.plot(Concentration_B_Data_1.iloc[:, 0], Concentration_B_Data_1.iloc[:, 1], color='red', label='% B')
ax2.set_ylabel('% Mobile Phase B', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(left=0, right=50)

# Save plot as PNG in the same directory as the first CSV file
output_dir = os.path.dirname(file_paths[0])
output_path = os.path.join(output_dir, 'Akta_Chromatogram_Comparison.png')
plt.savefig(output_path, bbox_inches='tight', dpi=300)
plt.close()