from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.encoding import data_dirs
import os

# Short script for training on an HPC
os.chdir("/xdisk/mtmarty/mtmarty/training/data")

eng = IsoDecEngine(phaseres=4)
topdirectory = "/xdisk/mtmarty/mtmarty/training/data"
dirs = [os.path.join(topdirectory, d) for d in data_dirs]
eng.create_merged_dataloader(dirs, "phase4", noise_percent=0, batchsize=32, double_percent=0)
eng.train_model(epochs=10)
eng.create_merged_dataloader(dirs, "phase4", noise_percent=0.0, batchsize=32, double_percent=0.2)
eng.train_model(epochs=10)
eng.create_merged_dataloader(dirs, "phase4", noise_percent=0.05, batchsize=32, double_percent=0)
eng.train_model(epochs=10)

eng.create_merged_dataloader(dirs, "phase4", noise_percent=0.02, batchsize=32, double_percent=0.05)
eng.train_model(epochs=10)
eng.train_model(epochs=10, lossfn="weightedcrossentropy")
eng.train_model(epochs=10, lossfn="focal")
eng.train_model(epochs=10)


