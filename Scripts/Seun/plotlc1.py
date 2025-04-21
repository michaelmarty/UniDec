from LCImporter import *

import os
import numpy as np
import matplotlib.pyplot as plt

def plot_data_from_paths(paths):
    # Initialize the maximum absorbance value to zero
    max_absorbance = 0

    for path in paths:
        datafiles = parse_allfiles(path)

        for i, datafile in enumerate(datafiles):
            data = datafile.get_data()

            # Find the maximum absorbance value across all chromatograms
            max_absorbance = max(max_absorbance, np.amax(data[:, 1]))

    for path in paths:
        datafiles = parse_allfiles(path)

        for i, datafile in enumerate(datafiles):
            info = datafile.get_info()
            data = datafile.get_data()

            # Normalize the y-axis to the maximum absorbance across all chromatograms
            #normalized_data = data[:, 1] / max_absorbance * 100


            # Create a separate plot for each chromatogram
            plt.plot(data[:, 0], data[:,1], label=f'Chromatogram {i + 1} - Absorbance')

            # Set axis labels and show legend for each plot
            plt.xlabel('Retention time (min)')
            #plt.ylabel('Normalized Absorbance')
            plt.ylabel("Absorbance (mAU)")
            plt.legend()

            # Set y-axis limits to 0 and 100
            plt.ylim(0, 100)
            plt.xlim(0,30)

            # Show the plot for the current chromatogram
            plt.title(f'LC Data - {info} - Chromatogram {i + 1}')

            # Display the plot instead of saving as a PDF file
            plt.show()

            # Close the current figure to clear the plot for the next iteration
            plt.close()

# Example usage
desktop_folder = os.path.join(os.path.expanduser("~"), "Desktop", "BCP_folder")
paths_to_plot = [
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\Blank.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\GFP.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\Blue Dextran.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-21 14-00-49\\GFP.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-14 16-13-52\\mRuby3.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_vary vol inject 2024-03-02 09-42-57\\010-D1F-A2-Blank.D",
    "Z:\\Group Share\\Seun\\LX-MS\LC\\LC Data\\OF_vary vol inject 2024-03-02 09-42-57\\026-D1F-E1-Blue Food dye.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-09 12-08-34\\20240309_Blank 10 ul.D",
    "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-09 12-08-34\\2024309_sf GFP 10 ul.D"

 # Add more paths as needed
]

# Create the output folder if it doesn't exist
os.makedirs(desktop_folder, exist_ok=True)

plot_data_from_paths(paths_to_plot)