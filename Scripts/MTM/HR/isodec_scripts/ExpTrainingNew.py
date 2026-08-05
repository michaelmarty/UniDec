from IsoDec.engine import *
import os
import time

# os.environ["CUDA_LAUNCH_BLOCKING"]="1"

startime = time.perf_counter()
# Define model

os.chdir("C:\\Data\\IsoNN\\exp_centroids")

eng = IsoDecEngine()

training_file = "exp_training_data_large.pth"
test_file = "exp_test_data_large.pth"

eng.create_training_dataloader(training_file, test_file)
eng.train_model(epochs=1, save=False)

endtime = time.perf_counter()
print("Time:", endtime - startime)
