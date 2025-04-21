import os
import glob
import matplotlib.pyplot as plt
from matplotlib import image as mpimg

# Specify the top directory containing the PNG files
topdir = "Z:\\Group Share\\Dash\\Ron - Orbi Native Data\\20241107 TD NAnodisc\\DPPG\\UniDec_Figures_and_Files"

# Get a list of all PNG files in the directory
png_files = glob.glob(os.path.join(topdir, '*.png'))

# Loop through each PNG file
for file_path in png_files:
    # Extract the file name without the extension
    file_name = os.path.basename(file_path).split('.')[0]

    # Read the image
    img = mpimg.imread(file_path)

    # Create a new figure
    plt.figure()

    # Display the image
    plt.imshow(img)

    # Add a title aligned top and center
    plt.title(file_name, loc='center', pad=20)

    # Remove axis
    plt.axis('off')

    # Save the updated figure
    plt.savefig(file_path, bbox_inches='tight', pad_inches=0.1)

    # Close the figure
    plt.close()