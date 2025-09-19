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


#def sigmoid(x, a, b, c, d):
    #y = d / (1.0 + np.exp(-a * (x - b))) + c
    #return y
def cube(x, a, b, c):
    y = a * x ** 2 + b * x + c  ##+ d * x**3 + e * x ** 4 + f * x ** 5
    return y

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"

xdata = np.array([0.273571103, 0.427924, 0.521518472, 0.606904818, 0.721550878, 0.818676778, 0.169139,  # batch a, 20
                  0.264946629, 0.377567405, 0.499840175, 0.611320635, 0.691793579, 0.791046772, 0.135019464,  # batch b
                  0.257864562, 0.402562646, 0.545817, 0.620564076, 0.662262, 0.789596279, 0.18042874])  # batch c, 40, 80
ydata = np.array([0.014557637, 0.039132, 0.076345411, 0.093313905, 0.106628672, 0.076113749, -0.02814,  # batch a, 20
                  0.009392632, 0.075266738, 0.115078044, 0.134266394, 0.104007865, 0.128838869, -0.058457228,  # batch b
                  0.034923576, 0.048626812, 0.114883, 0.138451348, 0.062596, 0.09148421, -0.024295153])  # batch c, 40, 80
ydata *= -1

popt, pcov = curve_fit(cube, xdata, ydata)
print(popt)

x = np.linspace(0, 1, 50)
y = cube(x, *popt)

rcParams['font.size'] = 16
#rcParams['text.usetex'] = True


ax = xdata[:7]
bx = xdata[7:14]
cx = xdata[14:]
ay = ydata[:7]
by = ydata[7:14]
cy = ydata[14:]

xav = np.mean([ax, bx, cx], axis=0)*100
xerr = np.std([ax, bx, cx], axis=0)*100
yav = np.mean([ay, by, cy], axis=0)*100
yerr = np.std([ay, by, cy], axis=0)*100

print(xav)
print(yav)

if False:
    for i, x in enumerate(xav):
        v1 = int(x)
        v2 = int(x + yav[i])
        v1, v2 = v2, v1
        print(v1, v2)
        savepath = os.path.join(dir, str(v1) + "_" + str(v2) + ".png")
        make_pie_chart(v1/100., v2/100., savefile=savepath)
    exit()


#plt.plot(xdata[:7]*100, ydata[:7] * 100, 'o', color="navy", label='Rep. 1')
#plt.plot(xdata[7:14] * 100, ydata[7:14] * 100, 'o', color="royalblue", label='Rep. 2')
#plt.plot(xdata[14:]*100, ydata[14:] * 100, 'o', color="steelblue", label='Rep. 3')

plt.plot(x * 100, y * 100, color="b", label='Fit', linestyle="-")

plt.plot(x*100, y*0, color="k", linestyle="--")
plt.errorbar(xav, yav, yerr=yerr, xerr=xerr, color="navy", linestyle=" ", capsize = 3)
plt.plot(xav, yav, color="yellow", marker="o", linestyle=" ", markeredgecolor='k')


plt.ylim(-20, 10)
plt.xlim(0,100)
plt.xlabel('X  (Total % POPC)')
plt.ylabel('X  \u2212  X  (\u0394 % POPC (AmtB \u2212 Total))')


ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#plt.legend(loc='best')
plt.tight_layout()


savepath = os.path.join(dir, "GlobalFit3_HR")
if False:
    plt.savefig(savepath+".png", transparent=True, dpi=1200)
    plt.savefig(savepath+".pdf", transparent=True, dpi=1200)
plt.show()
