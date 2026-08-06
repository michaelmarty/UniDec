from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.encoding_profile import data_dirs
import os
import time

def run_training(dirs, phaseres, centroid=True, phaselock=0, force_new=True, epochs=10):
    if centroid:
        centroidkey = "centroids"
        short='c_'
    else:
        centroidkey = "profiles"
        short='p_'

    ematkey="emat_" + short + str(phaseres) + "_" + str(phaselock)

    eng = IsoDecEngine(phaseres=phaseres)
    eng.create_merged_dataloader(dirs, f"profile", noise_percent=0, batchsize=32, double_percent=0.4,
                                 onedrop_percent=0.0, harmonic_percent=0, equilize=True, ematkey=ematkey,
                                 centroidkey=centroidkey)
    acc = eng.train_model(epochs=epochs, forcenew=force_new)
    return acc



topdirectory = "/groups/mtmarty/data"
topdirectory = "Z:\\Group Share\\JGP\\"

dirs = [os.path.join(topdirectory, d) for d in data_dirs]
dirs = ["Z:\\Group Share\\JGP\\PXD069439"]

# Short script for training on an HPC
os.chdir(topdirectory)

types = [[1, 0], [4, 0], [4, 1], [8, 0], [8, 1]]
epochs = 10
start_time = time.time()

accs = []

for t in types:
    phaseres = t[0]
    phaselock = t[1]
    print("Training: PhaseRes:", phaseres, "PhaseLock:", phaselock, "Centroid: True")
    t1 = time.perf_counter()
    acc1 = run_training(dirs, phaseres=phaseres, centroid=True, phaselock=phaselock, force_new=True, epochs=epochs)
    t1 = time.perf_counter() - t1
    print("Training: PhaseRes:", phaseres, "PhaseLock:", phaselock, "Centroid: False")
    t2 = time.perf_counter()
    acc2 = run_training(dirs, phaseres=phaseres, centroid=False, phaselock=phaselock, force_new=True, epochs=epochs)
    t2 = time.perf_counter() - t2

    accs.append([phaseres, phaselock, acc1*100, acc2*100, t1, t2])

# Print summary of accuracy results
print(f"For {epochs} epochs, the accuracy results are:")
headers = ["PhaseRes", "PhaseLock", "Centroid Accuracy", "Profile Accuracy", "Time Centroid (s)", "Time Profile (s)"]
row_fmt = "{:<10} {:<10} {:<20} {:<18} {:<18} {:<16}"
print(row_fmt.format(*headers))
print("-" * 100)
for phaseres, phaselock, acc1, acc2, t1, t2 in accs:
    print(row_fmt.format(phaseres, phaselock, f"{acc1:.4f}", f"{acc2:.4f}", f"{t1:.2f}", f"{t2:.2f}"))

print(f"Total time for training: {(time.time() - start_time)/60.:.2f} min")

# Save these result to a text file, make a tab separated table of the accuracy results
with open("accuracy_results.txt", "w") as f:
    f.write("\t".join(headers) + "\n")
    for phaseres, phaselock, acc1, acc2, t1, t2 in accs:
        f.write(f"{phaseres}\t{phaselock}\t{acc1:.4f}\t{acc2:.4f}\t{t1:.2f}\t{t2:.2f}\n")
