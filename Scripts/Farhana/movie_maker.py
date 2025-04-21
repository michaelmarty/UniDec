import glob
from PIL import Image
import os
import natsort
import numpy as np

def make_gif(frame_folder):
    os.chdir(frame_folder)
    files = glob.glob(f"{frame_folder}/*.bmp")
    files = natsort.natsorted(files)
    #files = files[100:250:3] + files[250:450] + files[450::3]
    frames = [Image.open(image) for image in files][::3]


    size = frames[0].size
    newsize = np.round(np.array(size) * 1).astype(int)
    newsize = tuple(newsize)
    print(size, newsize)
    frames = [f.resize(newsize) for f in frames]

    frame_one = frames[0]
    frame_one.save("output.gif", format="GIF", append_images=frames,
                   save_all=True, duration=85, loop=0)


if __name__ == "__main__":
    print("Starting")
    make_gif("Z:\\Group Share\\FAM\\VMD\\VMD\\Movie")
    print("Done")