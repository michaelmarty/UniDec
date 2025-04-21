import os
import glob
from pptx import Presentation
from pptx.util import Inches

# Specify the top directory containing the PNG files
topdir = "Z:\\Group Share\\Dash\\Ron - Orbi Native Data\\20241107 TD NAnodisc\\DPPC\\UniDec_Figures_and_Files"

# Get a list of all PNG files in the directory that end with "MassDefectFigure6"
png_files = glob.glob(os.path.join(topdir, '*MassDefectFigure6.png'))

# Create a PowerPoint presentation object
prs = Presentation()

# Loop through each PNG file
for file_path in png_files:
    # Add a slide with a blank layout
    slide_layout = prs.slide_layouts[6]  # 6 is for a blank slide
    slide = prs.slides.add_slide(slide_layout)

    # Add the image to the slide
    left = Inches(1)
    top = Inches(1)
    pic = slide.shapes.add_picture(file_path, left, top, width=Inches(8.5), height=Inches(6))

# Save the PowerPoint presentation
output_path = os.path.join(topdir, 'MassDefectFiguresPresentation.pptx')
prs.save(output_path)

print(f"PowerPoint presentation saved to {output_path}")