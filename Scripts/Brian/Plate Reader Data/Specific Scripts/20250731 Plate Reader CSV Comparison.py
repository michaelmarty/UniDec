# Import required libraries
import pandas as pd
import matplotlib.pyplot as plt
import os

# Set global font properties for all plot text
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# Define the path to the CSV data file
file_path = r'Z:\Group Share\BHT\Plate Reader Data\20250731\20250731_BHT_SF-GFP-MSP_postsec_ammoniumacetate_pH7.42_2M Sanitized.csv'

# Read the CSV file into a DataFrame
df = pd.read_csv(file_path, skiprows=0)

# Extract sample names from the first column (excluding header)
sample_names = df.iloc[1:, 0].values

# Remove the first column (sample names) from the DataFrame
df = df.iloc[:, 1:]

# Extract wavelength values from the first row
wavelengths = df.iloc[0].values

# Extract absorbance data for each sample (excluding header row)
absorbance_datasets = [row.values for _, row in df.iloc[1:].iterrows()]

# Get the directory of the data file for saving plots
data_dir = os.path.dirname(file_path)

# Indices for samples to compare (1st, 2nd, 5th, 6th, 9th, and 10th; 0-based indexing)
indices = [0, 1, 4, 5, 8, 9]

# Create the full-range plot
plt.figure(figsize=(16, 9))
for idx in indices:
    plt.plot(wavelengths, absorbance_datasets[idx], label=sample_names[idx])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Optical Density')
plt.legend()
plt.savefig(os.path.join(data_dir, "comparison_1_2_5_6_9_10.png"))
plt.close()

# Create the zoomed-in plot (300–500 nm, 0–0.5 OD)
plt.figure(figsize=(16, 9))
for idx in indices:
    plt.plot(wavelengths, absorbance_datasets[idx], label=sample_names[idx])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Optical Density')
plt.legend()
plt.xlim(300, 500)
plt.ylim(0, 0.5)
plt.savefig(os.path.join(data_dir, "comparison_1_2_5_6_9_10_zoom_300_500.png"))
plt.close()