import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Set plot font and size
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Load CSV data
file_path = r'Z:\Group Share\BHT\Plate Reader Data\20250731\20250731_BHT_SF-GFP-MSP_postsec_ammoniumacetate_pH7.42_2M Sanitized.csv'
df = pd.read_csv(file_path, skiprows=0)

# Extract sample names, wavelengths, and absorbance data
sample_names = df.iloc[1:, 0].values
wavelengths = df.iloc[0, 1:].values.astype(float)
absorbance_data = df.iloc[1:, 1:].astype(float).values

# Find index for 350 nm
idx_350 = (wavelengths == 350).nonzero()[0][0]
od_350 = absorbance_data[:, idx_350]

# Map sample indices to groups
group_map = {
    'Water Blank': [0],
    'Water': [1, 2, 3],
    '2M Ammonium Acetate Blank': [4],
    '2M Ammonium Acetate': [5, 6, 7],
    '4M Ammonium Sulfate Blank': [8],
    '4M Ammonium Sulfate': [9, 10, 11]
}

# Calculate mean and std for each group, set std to np.nan for blanks
group_names = []
group_means = []
group_stds = []
for group, indices in group_map.items():
    values = od_350[indices]
    group_names.append(group)
    group_means.append(np.mean(values))
    if 'Blank' in group:
        group_stds.append(np.nan)
    else:
        group_stds.append(np.std(values))

# Plot bar chart with error bars (no error bars for blanks)
plt.figure(figsize=(16, 12))  # Increased height
plt.bar(group_names, group_means, yerr=group_stds, capsize=8)
plt.xlabel('Sample Group')
plt.ylabel('Optical Density at 350 nm')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save figure
data_dir = os.path.dirname(file_path)
plt.savefig(os.path.join(data_dir, "bar_chart_OD_350nm_grouped_errorbars.png"))
plt.close()