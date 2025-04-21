import os
from pdf2image import convert_from_path

def pdf_to_png(pdf_path, output_folder):
    # Convert PDF to a list of images
    images = convert_from_path(pdf_path)

    # Save each image as a PNG file
    for i, image in enumerate(images):
        output_path = os.path.join(output_folder, f'{os.path.splitext(os.path.basename(pdf_path))[0]}_page_{i + 1}.png')
        image.save(output_path, 'PNG')
        print(f'Saved {output_path}')

def convert_all_pdfs_in_directory(topdir, output_folder):
    for root, _, files in os.walk(topdir):
        for file in files:
            if file.lower().endswith('.pdf'):
                pdf_path = os.path.join(root, file)
                pdf_to_png(pdf_path, output_folder)

if __name__ == "__main__":
    topdir = "Z:\\Group Share\\Dash\\Synapt Runs 20241218"  # Replace with your top directory path
    output_folder = "Z:\\Group Share\\Dash\\Synapt Runs 20241218"  # Replace with your output folder path
    os.makedirs(output_folder, exist_ok=True)
    convert_all_pdfs_in_directory(topdir, output_folder)