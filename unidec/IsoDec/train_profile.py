from pathlib import Path
import sys

path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.encoding_profile import data_dirs
import os
import time


def get_model_input_settings(phaseres, phaselock, mode):
    """Return the data keys and model filename for a profile-model input set."""
    if mode == "centroid":
        centroidkey = "centroids"
        short = "c_"
    elif mode == "profile":
        centroidkey = "profiles"
        short = "p_"
    elif mode == "both":
        centroidkey = "both"
        short = "b_"
    else:
        raise ValueError("Mode must be 'centroid', 'profile', or 'both'.")

    phaseres = int(phaseres)
    phaselock = int(phaselock)
    if phaseres not in (1, 4, 8):
        raise ValueError("Phase resolution must be 1, 4, or 8.")
    if phaselock not in (0, 1):
        raise ValueError("Phase lock must be 0 or 1.")
    ematkey = f"emat_{short}{phaseres}_{phaselock}"
    model_name = f"phase_model_{short}{phaseres}_{phaselock}.pth"
    return ematkey, centroidkey, model_name


def run_training(dirs, phaseres, mode="centroid", phaselock=0, force_new=True, epochs=10,
                 batchsize=32, double_percent=0.4):
    ematkey, centroidkey, _ = get_model_input_settings(phaseres, phaselock, mode)

    eng = IsoDecEngine(phaseres=phaseres)
    eng.create_merged_dataloader(dirs, f"profile", noise_percent=0, batchsize=batchsize,
                                 double_percent=double_percent,
                                 onedrop_percent=0.0, harmonic_percent=0, equilize=True, ematkey=ematkey,
                                 centroidkey=centroidkey)
    acc = eng.train_model(epochs=epochs, forcenew=force_new)
    return acc


if __name__ == "__main__":
    # TODO: new ROCs, timing check

    topdirectory = "/groups/mtmarty/data"
    topdirectory = "Z:\\Group Share\\JGP\\"

    dirs = [os.path.join(topdirectory, d) for d in data_dirs]
    # dirs = ["Z:\\Group Share\\JGP\\PXD069439"]

    # Short script for training on an HPC
    os.chdir(topdirectory)

    types = [[1, 0], [4, 0], [4, 1], [8, 0], [8, 1]]
    # types = [[1, 0]]
    modes = ["centroid", "profile"]
    modes = ["both"]
    epochs = 10
    double_percent = 0.4
    batch_size = 32
    start_time = time.time()

    accs = []

    for t in types:
        phaseres = t[0]
        phaselock = t[1]
        for m in modes:
            print(f"\n\nTraining: PhaseRes: {phaseres}, PhaseLock: {phaselock}, Mode: {m}")
            t1 = time.perf_counter()
            acc1 = run_training(dirs, phaseres=phaseres, mode=m, phaselock=phaselock, force_new=True, epochs=epochs,
                                double_percent=double_percent, batchsize=batch_size)
            t1 = time.perf_counter() - t1
            accs.append([phaseres, phaselock, m, acc1 * 100, t1])

    # Print summary of accuracy results
    print(f"\n\nFor {epochs} epochs, the accuracy results are:")
    headers = ["PhaseRes", "PhaseLock", "Mode", "Accuracy", "Time (s)"]
    row_fmt = "{:<10} {:<10} {:<10} {:<18} {:<18}"
    print(row_fmt.format(*headers))
    print("-" * 100)
    for phaseres, phaselock, m, acc1, t1 in accs:
        print(row_fmt.format(phaseres, phaselock, m, f"{acc1:.4f}", f"{t1:.2f}"))

    print(f"Total time for training: {(time.time() - start_time) / 60.:.2f} min")

    os.chdir(topdirectory)
    # Save these result to a text file, make a tab separated table of the accuracy results
    with open("accuracy_results_both.txt", "w") as f:
        f.write("\t".join(headers) + "\n")
        for phaseres, phaselock, m, acc1, t1 in accs:
            f.write(f"{phaseres}\t{phaselock}\t{m}\t{acc1:.4f}\t{t1:.2f}\n")
