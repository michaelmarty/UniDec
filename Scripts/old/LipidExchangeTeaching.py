import pymzml
import os
import numpy as np
import matplotlib.pyplot as plt

file = "MTM150515_PSPC_50.mzML"
dir = "C:\\Users\\Michael Marty\\Box Sync\\Teaching\\525 Fall 2016\\Data"
path = os.path.join(dir, file)
msrun = pymzml.run.Reader(path)

times = []
tic = []

for spectrum in msrun:
    try:
        time = float(spectrum["scan start time"])
    except:
        break
    mz = np.array(spectrum.mz)
    intensity = np.array(spectrum.i)
    times.append(time)
    tic.append(np.sum(intensity))
    print("Retention Time: ", time, "Peak m/z value: ", mz[np.argmax(intensity)], "Maximum Intensity: ", np.amax(
        intensity))

plt.plot(times, tic)
plt.ylabel("Total Ion Current")
plt.xlabel("Retention Time")
plt.show()
