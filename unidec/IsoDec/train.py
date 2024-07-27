from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
import os

# Short script for training on an HPC
os.chdir("/xdisk/mtmarty/mtmarty/training/data")

eng = IsoDecEngine()
topdirectory = "/xdisk/mtmarty/mtmarty/training/data"
dirs = [os.path.join(topdirectory, d) for d in data_dirs]
eng.create_merged_dataloader(dirs, "phase82", noise_percent=0.2, batchsize=32, double_percent=0.2)

eng.train_model(epochs=60)
