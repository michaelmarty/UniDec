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

xdata = np.array([0.140848709, 0.288593194, 0.388878424, 0.478126487, 0.584394742, 0.690072965, 0.798087733,  # batch a
                  0.176108533, 0.296031736, 0.393540509, 0.502495335, 0.565091666, 0.70119642, 0.812863099,  # batch b
                  0.154151584, 0.282013477, 0.350053079, 0.493714773, 0.599015074, 0.705821282, 0.783380951])  # batch c
ydata = np.array([-0.001449824, -0.00864566, 0.004648944, -0.022173523, 0.01761642, 0.016373769, 0.02639904,  # batch a
                  -0.032142842, -0.004353371, 0.016886864, - 0.033814336, 0.023037043, 0.019856835, 0.01384501,  # batch b
                  -0.038826509, -0.042101578, 0.012053613, 0.011109559, - 0.014327672, 0.002162886, 0.02379928])  # batch c
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


savepath = os.path.join(dir, "GlobalFit3_HR_Control2")
if False:
    plt.savefig(savepath+".png", transparent=True, dpi=1200)
    plt.savefig(savepath+".pdf", transparent=True, dpi=1200)
plt.show()
