from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
import os


os.chdir("/xdisk/mtmarty/mtmarty/training")
eng = IsoDecEngine(1)
eng.create_training_dataloader("large32x2")
#eng.create_training_dataloader("exp_training_data_small.pth", "exp_test_data_small.pth")
eng.train_model(epochs=120)
