import glob
from PIL import Image
import os
import natsort


def make_gif(frame_folder):
    files = glob.glob(os.path.join(frame_folder, "*.png"))
    files = natsort.natsorted(files)

    frames = [Image.open(image) for image in files]

    if not frames:
        print("No PNG files found.")
        return

    frames[0].save(
        os.path.join(frame_folder, "my_awesome.gif"),
        format="GIF",
        append_images=frames[1:],  # prevent duplicating first frame
        save_all=True,
        duration=45,
        loop=0
    )


if __name__ == "__main__":
    make_gif(r"C:\Python\Animations")