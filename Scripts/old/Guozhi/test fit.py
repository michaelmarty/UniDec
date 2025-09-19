import numpy as np
import pylab
from scipy.optimize import curve_fit
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import os

def make_pie_chart(inner, outer, savefile=None):
    fig = plt.figure(figsize=(1.0,1.0))
    ax = plt.gca()
    ax.set_position(pos=[0,0,1,1])
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
    #plt.tight_layout(h_pad=0.01, w_pad=0.01)
    #plt.show()


def sigmoid(x, a, b, c, d):
    y = d / (1.0 + np.exp(-a * (x - b))) + c
    return y


dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"

xdata = np.array([0.273571103, 0.383556608, 0.521518472, 0.606904818, 0.721550878, 0.818676778, 0.173656229,  # batch a
                  0.264946629, 0.377567405, 0.499840175, 0.611320635, 0.691793579, 0.791046772, 0.135019464,  # batch b
                  0.257864562, 0.402562646, 0.489138882, 0.620564076, 0.703541189, 0.789596279, 0.18042874,  # batch c
                  0.265460765, 0.387895553, 0.503499176, 0.612929843, 0.705628549, 0.799773276, 0.143942109])  # batch a repeat
ydata = np.array([0.014557637, 0.127916528, 0.076345411, 0.093313905, 0.106628672, 0.076113749, -0.012603438,   # batch a
                  0.009392632, 0.075266738, 0.115078044, 0.134266394, 0.104007865, 0.128838869, -0.058457228,  # batch b
                  0.034923576, 0.048626812, 0.184739496, 0.138451348, 0.182809889, 0.091484210, -0.024295153, # batch c
                  0.019624615, 0.083936693, 0.125387650, 0.122010549, 0.131148809, 0.098812276, -0.011494427])  # batch a repeat

ydata *= -1

popt, pcov = curve_fit(sigmoid, xdata, ydata)
print(popt)

x = np.linspace(0, 1, 50)
y = sigmoid(x, *popt)

rcParams['font.size'] = 16
#rcParams['text.usetex'] = True


ax = xdata[:7]
bx = xdata[7:14]
cx = xdata[14:21]
dx = xdata[21:]

ay = ydata[:7]
by = ydata[7:14]
cy = ydata[14:21]
dy = xdata[21:]

xav = np.mean([ax, bx, cx, dx], axis=0)*100
xerr = np.std([ax, bx, cx, dx], axis=0)*100
yav = np.mean([ay, by, cy, dy], axis=0)*100
yerr = np.std([ay, by, cy, dy], axis=0)*100

print(xav)
print(yav)

plt.plot(x * 100, y * 100, color="b", label='Fit', linestyle="-")

plt.plot(x*100, y*0, color="k", linestyle="--")
plt.errorbar(xav, yav, yerr=yerr, xerr=xerr, color="navy", linestyle=" ", capsize = 3)
plt.plot(xav, yav, color="yellow", marker="o", linestyle=" ", markeredgecolor='k')


plt.ylim(-20, 10)
plt.xlim(0,100)
plt.xlabel('Total % POPC')
plt.ylabel('\u0394 % POPC (AmtB \u2212 Total)')