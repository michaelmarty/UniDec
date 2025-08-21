import pandas as pd
import matplotlib.pyplot as plt
import os  # Import os for file path operations

# Define the path to the CSV file
file_path = r'Z:\Group Share\BHT\Plate Reader Data\20250731\20250731_BHT_SF-GFP-MSP_postsec_ammoniumacetate_pH7.42_2M Sanitized.csv'

# Read the CSV file into a DataFrame
df = pd.read_csv(file_path, skiprows=0)

# Extract sample names from the first column, excluding the header row
sample_names = df.iloc[1:, 0].values

# Remove the first column (sample names) for plotting purposes
df = df.iloc[:, 1:]

# Get the wavelength values from the first row
wavelengths = df.iloc[0].values

# Get absorbance data for each sample (rows after the first)
absorbance_datasets = [row.values for _, row in df.iloc[1:].iterrows()]

# Get the directory of the data file to save plots
data_dir = os.path.dirname(file_path)

# Loop through each sample and save its plot as a PNG in the data directory
for i, absorbance in enumerate(absorbance_datasets):
    plt.figure()
    plt.plot(wavelengths, absorbance)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Optical Density')
    plt.title(f'Sample: {sample_names[i]}')
    plt.savefig(os.path.join(data_dir, f"{sample_names[i]}.png"))  # Save plot
    plt.close()  # Close the figure to free memory