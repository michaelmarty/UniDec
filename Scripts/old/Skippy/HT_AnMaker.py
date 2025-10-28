import glob
from PIL import Image
import os
import natsort


def make_gif(frame_folder):
    os.chdir(frame_folder)
    files = glob.glob(f"{frame_folder}/*.png")
    files = natsort.natsorted(files)
    files = files[100:250:3] + files[250:450] + files[450::3]
    frames = [Image.open(image) for image in files]
    frame_one = frames[0]
    frame_one.save("my_awesome.gif", format="GIF", append_images=frames,
                   save_all=True, duration=45, loop=0)


if __name__ == "__main__":
    make_gif("C:\\Python\\Animations\\")
