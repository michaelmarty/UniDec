import numpy as np
import matplotlib.pyplot as plt
import os

file1 = "D:\\UArizona\\Game On\\TA\\Fall 2024\\testtubes\\testtubes\\System 1_ td mix 300ppm rep03.ASC"
ydata = np.loadtxt(file1, skiprows=25, usecols=0, delimiter=",", comments = "I")
ydata1 = ydata / 1000

rate = 5.0
xdata1 = np.arange (0, len(ydata1))/float(rate)/60


def datachop(x, y, minx, maxx):
    boo1 = np.logical_and( x <= maxx, x >= minx)
    return x[boo1], y[boo1]

# Define three different minx and maxx pairs
ranges = [
    (1.21, 1.70),  # Peak 1
    (1.709, 2.6),  # Peak 2
    (4.1, 5.3)   # Peak 3
]

areas = []  # To store areas of each peak

# Plotting setup
plt.figure(figsize=(10, 6))

# Plot the full dataset
plt.plot(xdata1, ydata1, label='Full Data', color='black', linewidth=1.5)

for minx, maxx in ranges:
    cutx, cuty = datachop(xdata1, ydata1, minx, maxx)
    area = np.trapz(cuty, x=cutx)
    areas.append(area)

    # Plot the chopped data
    plt.fill_between(cutx, cuty, alpha=0.3)
    plt.plot(cutx, cuty, label=f'Area = {area:.2f} from {minx} to {maxx}', linewidth =2)


# Final plot settings
plt.xlabel("Time (min)")
plt.ylabel("Intensity (mV)")
plt.title("HPLC Chromatogram with Integrated Areas")
plt.legend()
plt.show()


