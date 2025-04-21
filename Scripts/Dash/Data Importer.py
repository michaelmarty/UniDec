import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
import glob
from matplotlib import rcParams
import numpy as np


rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 18
rcParams['font.weight'] = 700
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.size'] = 4
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.size'] = 4
rcParams['ytick.major.width'] = 2
rcParams['figure.facecolor'] = 'white'
rcParams['figure.edgecolor'] = 'white'
rcParams['figure.autolayout'] = True

#plt.ion()




# Define the directory containing the raw files
topdir = "Z:\\Group Share\\Dash\\Orbi_LCMS_Test"
# Get a list of all .raw files in the directory
files = glob.glob(f'{topdir}/*.raw', recursive=False)

# Print the found files
print(files)

# Prepare a figure for the plot
plt.figure(figsize=(16, 9))

# Initialize a list to hold all data for batch plotting
all_data = []

# Loop through each file and extract the TIC data
for path in files:
    print(f"Processing file: {path}")
    try:
        # Create the importer for the current file
        importer = ImporterFactory.create_importer(path)

        # Extract the Total Ion Chromatogram (TIC) data
        data = importer.get_tic()

        # Collect data for batch processing
        all_data.append((data[:, 0], data[:, 1], os.path.basename(path)))

    except Exception as e:
        print("File Failed:", path)
        print(e)

# Define the range for filling the area under the curves
fill_start = 10.3
fill_end = 10.8

# Plot all data at once and fill the area under the curve
for time, intensity, label in all_data:
    plt.plot(time, intensity, label=label)

    # Fill the area under the curve between the specified limits
    plt.fill_between(time, intensity, where=((time >= fill_start) & (time <= fill_end)),
                     alpha=0.3)  # Adjust alpha for transparency


plt.xlabel('Time (min)')
plt.ylabel('Intensity')
plt.title('Superimposed TIC from Multiple Files with Area Fill')
plt.tight_layout()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)


# Show the plot
plt.show()