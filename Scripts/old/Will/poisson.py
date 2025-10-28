import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson as p

mu=1

xvals=np.arange(0,mu*4)

yvals=np.array([p.pmf(x,mu) for x in xvals])

plt.bar(xvals, yvals)
plt.show()