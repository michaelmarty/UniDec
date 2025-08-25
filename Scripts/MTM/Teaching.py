import os
import numpy as np
import matplotlib.pyplot as plt

dir = r"C:\Data\EIC Dataset"
os.chdir(dir)

files = os.listdir(dir)

for file in files:
    if file.endswith("example.txt"):
        data = np.loadtxt(file)
        plt.plot(data[:, 0], data[:, 1], label=file)

plt.legend()
plt.show()
