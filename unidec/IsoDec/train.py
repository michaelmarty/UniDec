from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.encoding import data_dirs
import os


topdirectory = "/groups/mtmarty/data"
dirs = [os.path.join(topdirectory, d) for d in data_dirs]

# Short script for training on an HPC
os.chdir(topdirectory)

eng = IsoDecEngine(phaseres=8)
eng.create_merged_dataloader(dirs, "phase83", noise_percent=0, batchsize=32, double_percent=0.4,
                             onedrop_percent=0.0, harmonic_percent=0, equilize=True)
eng.train_model(epochs=10, forcenew=True)

eng = IsoDecEngine(phaseres=1)
eng.create_merged_dataloader(dirs, "phase1", noise_percent=0, batchsize=32, double_percent=0.4,
                             onedrop_percent=0.0, harmonic_percent=0, equilize=True)
eng.train_model(epochs=10, forcenew=True)
#eng.create_merged_dataloader(dirs, datatype, noise_percent=0.0, batchsize=32, double_percent=0.2)
#eng.train_model(epochs=10)
#eng.create_merged_dataloader(dirs, datatype, noise_percent=0.05, batchsize=32, double_percent=0)
#eng.train_model(epochs=10)

#eng.create_merged_dataloader(dirs, datatype, noise_percent=0.02, batchsize=32, double_percent=0.05)
#eng.train_model(epochs=10)
#eng.train_model(epochs=10, lossfn="weightedcrossentropy")
#eng.train_model(epochs=10, lossfn="focal")
#eng.train_model(epochs=10)
#eng.train_model(epochs=10)


