import numpy as np
import matplotlib.pyplot as plt
import os

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"

def simple_pie_chart(ratio, savefile=None):
    fig, ax = plt.subplots()
    inner_colors = ["darkslategray", "darkgreen"]

    ax.pie([ratio, 1 - ratio], radius=1, colors=inner_colors,
           wedgeprops=dict(width=1, edgecolor='w'))
    ax.set(aspect="equal", title='')
    if savefile is not None:
        plt.savefig(savefile, transparent=True)
    plt.show()

#savepath = os.path.join(dir, "single_27.png")
#simple_pie_chart(.27, savefile=savepath)
#exit()

def make_pie_chart(inner, outer, savefile=None):
    fig, ax = plt.subplots()
    size = 0.75

    outer_colors = ["slategray", "forestgreen"]
    inner_colors = ["darkslategray", "darkgreen"]

    ax.pie([outer, 1-outer], radius=1, colors=outer_colors,
           wedgeprops=dict(width=size, edgecolor='w'))

    ax.pie([inner, 1-inner], radius=size, colors=inner_colors,
           wedgeprops=dict(width=size, edgecolor='w'))

    ax.set(aspect="equal", title='')
    if savefile is not None:
        plt.savefig(savefile, transparent=True)
    plt.show()

v1=24
v2=v1+0
savepath = os.path.join(dir, str(v1)+"_"+str(v2)+".png")
print(savepath)
make_pie_chart(v1/100., v2/100, savefile=savepath)

exit()
path = "Z:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\crude titration"
os.chdir(path)

file1 = "total lipid ratio.txt"
file2 = "AmtB ND.txt"

totdat = np.loadtxt(file1)
amdat = np.loadtxt(file2)

amdat = np.transpose([amdat, 1-amdat])
totdat = np.transpose([totdat, 1-totdat])

i=0

for i in range(0, len(amdat), 1):
    fig, ax = plt.subplots()
    size = 0.75

    cmap = plt.get_cmap("tab20c")
    outer_colors = ["slategray", "forestgreen"]
    inner_colors = ["darkslategray", "darkgreen"]

    ax.pie(totdat[i], radius=1, colors=outer_colors,
           wedgeprops=dict(width=size, edgecolor='w'))

    ax.pie(amdat[i], radius=size, colors=inner_colors,
           wedgeprops=dict(width=size, edgecolor='w'))

    ax.set(aspect="equal", title='')
    dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
    path = os.path.join(dir, "figure_"+str(i)+".png")
    plt.savefig(path, transparent=True)
plt.show()