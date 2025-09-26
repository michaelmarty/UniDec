import numpy as np
import matplotlib.pyplot as plt


def Gaussian(xvalues, mean, standard_deviation):
    sdsquared = standard_deviation**2
    exponent = -(xvalues-mean)**2/(2*sdsquared)
    preexponential = 1 /np.sqrt(2 * np.pi * sdsquared)
    return preexponential * np.exp(exponent)
#defines a gaussian function

x=np.arange(500,1100)
# sets the x range
y= Gaussian(x,10000,20)
z= Gaussian(x,8000,20)
w= Gaussian(x,6666.667,20)
# mass at 4000kDa charge states 4,5,6
# intensities 2500, 5000, 1666.667

plt.plot(x,w)
plt.plot(x,z)
plt.plot(x,y)
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.show()
