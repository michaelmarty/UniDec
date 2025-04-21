import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Gradient Color Bar Plots
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import colors as mcolors, path

file1 ="Z:\Group Share\Dash\TOC input.txt"
input1 = np.loadtxt(file1)

x = input1[:,0]
y = input1[:,1]



def gradientbars(bars, ydata, cmap):
    ax = bars[0].axes
    lim = ax.get_xlim() + ax.get_ylim()
    ax.axis(lim)
    for bar in bars:
        bar.set_facecolor("none")
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        grad = np.atleast_2d(np.linspace(0, 1 * h / max(ydata), 256)).T
        ax.imshow(grad, extent=[x, x + w, y, y + h], origin='lower', aspect="auto",
                  norm=cm.colors.NoNorm(vmin=0, vmax=1), cmap=plt.get_cmap(cmap))


fig, ax = plt.subplots()
my_bar = ax.bar(x, y)
gradientbars(my_bar, y, 'viridis')

plt.savefig("Z:\Group Share\Dash\TOC Output.pdf")

plt.show()