import numpy as np
import pylab
from scipy.optimize import curve_fit

def sigmoid(x, a, b, c, d):
     y = d / (1.0 + np.exp(-a*(x-b)))+c
     return y

xdata = np.array([0.939780263, 0.836706233, 0.67277635, 0.545816804, 0.350862464, 0.246023175, 0.189819536, 0.132491474, 0.042712802])
ydata = np.array([0.072391, 0.087088, 0.10168, 0.114883, -0.02115, -0.02783, -0.04156, -0.10999, -0.09633])

popt, pcov = curve_fit(sigmoid, xdata, ydata)
print(popt)

x = np.linspace(0, 1, 50)
y = sigmoid(x, *popt)

pylab.plot(xdata, ydata, '^', color="navy", label='Data')
pylab.plot(x,y, label='Fit', color="b", linestyle="-")
pylab.ylim(-0.15, 0.15)
pylab.xlabel('Total Final % PC')
pylab.ylabel('% PC Difference')
pylab.legend(loc='best')
pylab.show()


