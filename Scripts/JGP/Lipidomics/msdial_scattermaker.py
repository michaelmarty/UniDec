import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from Scripts.LipiDec.AddFeaturestoDF import *
from matplotlib.lines import Line2D

# Set backend to Agg
mpl.use('WxAgg')

file = "Z:\\Group Share\\Annika\\Stellar\\EC\\Raw Data\\Rt_1_2025_09_22_11_28_09.csv"

classes_to_ignore = ["RIKEN"]

df = pd.read_csv(file, delimiter=",")

unique_classes = []
classes = []
rts = []
mzs = []

for i in range(len(df)):
    name = df["Metabolite name"][i]
    curr_class = get_class(name)
    if curr_class in classes_to_ignore:
        continue
    classes.append(curr_class)
    rts.append(df["Average Rt(min)"][i])
    mzs.append(df["Average Mz"][i])
    if curr_class not in unique_classes:
        unique_classes.append(curr_class)



colors = get_colors(unique_classes)
class_colors = {}
for i in range(len(unique_classes)):
    class_colors[unique_classes[i]] = colors[i]

custom_legend = []
for i in range(len(unique_classes)):
    custom_legend.append(Line2D([0], [0], marker='o', color=colors[i], label=unique_classes[i],
                                markerfacecolor=colors[i], linestyle='None', markersize=10))


for i in range(len(classes)):
    plt.scatter(rts[i], mzs[i], color=class_colors[classes[i]])
plt.legend(handles=custom_legend)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Retention Time (min)", fontsize=16)
plt.ylabel("m/z", fontsize=16)
plt.show()

