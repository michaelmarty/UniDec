import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
import os
from multipass_ims import run_fft

datadir = "C:\Data\Varun\BSA_z17"
os.chdir(datadir)

files = ud.match_files(datadir, "*.csv")

# Loop through the files, extract the data, and add it together into a single array
data = np.zeros((0, 2))
for file in files:
    skip = ud.header_test(file, delimiter=",")
    newdata = np.loadtxt(file, delimiter=",", skiprows=skip)
    try:
        max_x_data = np.amax(data[:, 0])
    except Exception:
        max_x_data = 0

    newxaxis = np.append(data[:, 0], newdata[:, 0][newdata[:, 0] > max_x_data])
    newyaxis = np.zeros_like(newxaxis)

    newyaxis2 = ud.mergedata(np.transpose([newxaxis, newyaxis]), newdata)
    newyaxis += newyaxis2[:, 1]
    newyaxis[:len(data[:, 1])] += data[:, 1]

    data = np.transpose([newxaxis, newyaxis])

    # data = np.append(data, newdata, axis=0)

fftdata = run_fft(data, mintime=20, maxtime=60)

# Print the maximum peak in the fftdata
maxindex = np.argmax(fftdata[:, 1])
maxtime = fftdata[maxindex, 0]
print(maxtime)

# Plot the array
plt.figure()
plt.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]))
plt.plot(fftdata[:, 0], fftdata[:, 1] / np.amax(fftdata[:, 1]))

plt.xlabel("Time (ms)")
plt.vlines(maxtime, 0, fftdata[maxindex, 1], linestyles='--', color="k")
plt.text(maxtime, fftdata[maxindex, 1], str(maxtime))
plt.show()
