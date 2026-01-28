import os
import pandas as pd
import matplotlib.pyplot as plt

# Set plot font properties for consistency
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# List of file paths to chromatogram data
file_paths = [
    r'Z:\Group Share\BHT\Akta FPLC Data\20250807\20250807_BHT_anionexchange_GroEL_35mLsample 001.csv',
    r'Z:\Group Share\BHT\Akta FPLC Data\20250807\20250807_BHT_anionexchange_GroEL_35mLsample 002.csv'
]
# Labels and colors for each dataset
labels = ['Run 1', 'Run 2']
colors = ['blue', 'green']

# Read CSV files, skipping first 2 rows, using UTF-16-LE encoding and tab separator
dfs = [pd.read_csv(fp, skiprows=2, encoding="utf-16-le", sep='\t') for fp in file_paths]
# Extract first two columns (assumed to be volume and UV absorbance)
UV_datas = [df.iloc[:, :2] for df in dfs]
# Extract columns 4 and 5 from the first file (assumed to be concentration data)
Concentration_B_Data_1 = dfs[0].iloc[:, 4:6]

# Create a figure and primary axis for UV data
fig, ax1 = plt.subplots(figsize=(16, 9))
for UV_Data, label, color in zip(UV_datas, labels, colors):
    # Plot UV absorbance vs. volume for each dataset
    ax1.plot(UV_Data.iloc[:, 0], UV_Data.iloc[:, 1], color=color, label=label)
ax1.set_xlabel('Volume (mL)')
ax1.set_ylabel('UV Absorbance (mAU)')
ax1.tick_params(axis='y')
ax1.set_xlim(left=0, right=320)
ax1.set_ylim(bottom=-5)
ax1.legend(loc='upper right', fontsize=18)

# Create a secondary y-axis for mobile phase B concentration
ax2 = ax1.twinx()
ax2.plot(Concentration_B_Data_1.iloc[:, 0], Concentration_B_Data_1.iloc[:, 1], color='red', label='% B')
ax2.set_ylabel('% Mobile Phase B', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(left=0, right=320)

# Save the figure to the same directory as the input files
output_dir = os.path.dirname(file_paths[0])
output_path = os.path.join(output_dir, 'Akta_Chromatogram_Comparison.png')
plt.savefig(output_path, bbox_inches='tight', dpi=300)
plt.close()