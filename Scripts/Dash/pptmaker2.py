import os
import glob
from pptx import Presentation
from pptx.util import Inches

# Specify the subdirectory containing the PNG files
subdir = "Z:\\Group Share\\Dash\\Ron - Orbi Native Data\\20241107 TD NAnodisc\\DPPC\\UniDec_Figures_and_Files\\SelectedFigures"

# Get a list of all PNG files in the subdirectory
png_files = glob.glob(os.path.join(subdir, '*.png'))


# Function to group files by their base name (omitting the last word separated by an underscore)
def group_files(files):
    grouped = {}
    for file_path in files:
        base_name = '_'.join(os.path.basename(file_path).split('_')[:-1])
        if base_name not in grouped:
            grouped[base_name] = []
        grouped[base_name].append(file_path)
    return grouped


# Group the PNG files
grouped_files = group_files(png_files)

# Create a PowerPoint presentation object
prs = Presentation()

# Loop through each group of files
for base_name, files in grouped_files.items():
    # Add slides with two images per slide
    for i in range(0, len(files), 2):
        slide_layout = prs.slide_layouts[6]  # 6 is for a blank slide
        slide = prs.slides.add_slide(slide_layout)

        # Add the first image
        left = Inches(1)
        top = Inches(1)
        pic = slide.shapes.add_picture(files[i], left, top, width=Inches(4.25), height=Inches(3))

        # Add the second image if it exists
        if i + 1 < len(files):
            left = Inches(5.5)
            top = Inches(1)
            pic = slide.shapes.add_picture(files[i + 1], left, top, width=Inches(4.25), height=Inches(3))

# Save the PowerPoint presentation
output_path = os.path.join(subdir, 'GroupedFiguresPresentation.pptx')
prs.save(output_path)

print(f"PowerPoint presentation saved to {output_path}")