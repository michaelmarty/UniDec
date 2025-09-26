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
#def cube(x, a, b, c):
    #y = a * x ** 2 + b * x + c  ##+ d * x**3 + e * x ** 4 + f * x ** 5
    #return y

def linear (x, m, b):
    y = m * x + b
    return y

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"

xdata = np.array([0.140848709, 0.214688221, 0.374726296, 0.451188634, 0.568101044, 0.675784073, 0.769187668,  # batch a
                  0.176108533, 0.220167732, 0.385176052, 0.50730746, 0.557696086, 0.696266512, 0.780475773,  # batch b
                  0.154151584, 0.245799469, 0.351611912, 0.502834049, 0.608601154, 0.717469657, 0.77586694])  # batch c
ydata = np.array([-0.001449824, 0.000769344, -0.004401329, 0.00445419, 0.009396475, 0.028437854, 0.02084548,  # batch a
                  -0.032142842, 0.001744203, 0.008522407, 0.036487872, 0.016685921, 0.027136493, 0.024174084,  # batch b
                  -0.038826509, -0.032524754, 0.062508231, -0.021740258, 0.023003149, 0.035930076, 0.013001233])  # batch c
ydata *= -1

print(np.mean(ydata)*100, np.std(ydata)*100)
n = np.round(np.abs(ydata)*100)
print(np.mean(n), np.std(n))

popt, pcov = curve_fit(linear, xdata, ydata)
print(popt)

x = np.linspace(0, 1, 50)
y = linear(x, *popt)

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

#plt.plot(x * 100, y * 100, color="b", label='Fit', linestyle="-")

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


savepath = os.path.join(dir, "GlobalFit3_HR_Control")
if True:
    plt.savefig(savepath+".png", transparent=True, dpi=1200)
    plt.savefig(savepath+".pdf", transparent=True, dpi=1200)
plt.show()
