import os
import glob
import shutil

# Specify the top directory containing the PNG files
topdir = "Z:\\Group Share\\Dash\\Ron - Orbi Native Data\\20241107 TD NAnodisc\\DPPG\\UniDec_Figures_and_Files"

# Specify the subdirectory name
subdir = os.path.join(topdir, 'SelectedFigures')

# Create the subdirectory if it doesn't exist
os.makedirs(subdir, exist_ok=True)

# List of patterns to match
patterns = ["*MassDefectFigure6.png", "*MassDefectFigure1.png", "*MassDefectFigure2.png",
            "*MetaFigure1.png", "*MetaFigure2.png", "*MetaFigure6.png"]

# Loop through each pattern and copy matching files to the subdirectory
for pattern in patterns:
    for file_path in glob.glob(os.path.join(topdir, pattern)):
        shutil.copy(file_path, subdir)

print(f"Selected PNG files have been copied to {subdir}")