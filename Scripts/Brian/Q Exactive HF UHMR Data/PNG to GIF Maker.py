import os
from PIL import Image

def create_gif_from_pngs(folder_path, output_gif, duration=500):
    """
    Create an animated GIF from a series of PNG files in a specified folder.

    Parameters:
        folder_path (str): Path to the folder containing PNG files.
        output_gif (str): Path to save the output GIF file.
        duration (int): Duration for each frame in milliseconds. Default is 500ms.
    """
    # Get and sort PNG files numerically if filenames contain numbers
    png_files = sorted(
        (f for f in os.listdir(folder_path) if f.endswith('.png')),
        key=lambda x: int(''.join([char for char in x if char.isdigit()])) if any(char.isdigit() for char in x) else x
    )

    if not png_files:
        print("No PNG files found in the specified folder.")
        return

    # Load images
    images = [Image.open(os.path.join(folder_path, file)) for file in png_files]

    # Save as GIF with only one frame displayed at a time
    images[0].save(
        output_gif,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0,
        disposal=2  # Ensures only one frame is displayed at a time
    )
    print(f"Animated GIF saved as {output_gif}")

if __name__ == "__main__":
    folder_path = input("Enter the folder path containing PNG files: ").strip('"')
    output_gif = os.path.join(folder_path, "output.gif")
    duration = int(input("Enter the duration for each frame in milliseconds (default 500): ") or 500)

    create_gif_from_pngs(folder_path, output_gif, duration)
