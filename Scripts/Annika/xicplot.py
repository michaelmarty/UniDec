import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LinearSegmentedColormap

directory = "Z:\\Group Share\\Annika\\Figures"
os.chdir(directory)

file1 = "HEK_pos_TIC.txt"
file2 = "HEK_neg_TIC.txt"

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14

# --- CUSTOM COLORMAPS (NEW) ---
warm_cmap = LinearSegmentedColormap.from_list(
    "pink_orange_red",
    ["#ff5fa2", "#ff9f45", "#d62828"]
)

cool_cmap = LinearSegmentedColormap.from_list(
    "purple_blue_green",
    ["#6a00f4", "#1e88e5", "#2ecc71"]
)


fig, ax = plt.subplots(figsize=(11, 5))

# --- GRADIENT LINE FUNCTION (NEW) ---
def gradient_line(x, y, cmap, linewidth=2):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = Normalize(x.min(), x.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(x)
    lc.set_linewidth(linewidth)
    return lc

# --- ESI+ (WARM GRADIENT) ---
lc1 = gradient_line(data1[:, 0], data1[:, 1], warm_cmap, linewidth=2)
ax.add_collection(lc1)

# --- ESI− (COOL GRADIENT) ---
lc2 = gradient_line(data2[:, 0], data2[:, 1], cool_cmap, linewidth=2)
ax.add_collection(lc2)

# --- REQUIRED FOR VISIBILITY (NEW) ---
ax.autoscale()
ax.margins(x=0)

# --- ORIGINAL AXIS STYLING (UNCHANGED) ---
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params(axis='both', width=2, length=5)

ax.set_xlim(0, 20)
ax.set_ylim(0, None)
ax.set_xlabel("Retention Time (min)", fontsize=18)
ax.set_ylabel("Intensity", fontsize=18)

ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax.set_xticks(ax.get_xticks())
ax.set_yticks(ax.get_yticks())
ax.tick_params(labelsize=18)

# --- MANUAL LEGEND (NEW, REQUIRED) ---
from matplotlib.lines import Line2D
legend_lines = [
    Line2D([0], [0], color=warm_cmap(0.8), lw=2, label="ESI+"),
    Line2D([0], [0], color=cool_cmap(0.8), lw=2, label="ESI−")
]
ax.legend(handles=legend_lines)

plt.tight_layout()
plt.show()
