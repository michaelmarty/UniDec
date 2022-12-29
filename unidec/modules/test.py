import numpy as np
import matplotlib.pyplot as plt
from unidec import tools as ud
from copy import deepcopy
import warnings

warnings.filterwarnings("error")

x = np.arange(0., 10., 0.1)

mu = 1
y = ud.ndis(x, 5, mu)
mu2 = 2
p = ud.ndis(x, 0, mu2) + ud.ndis(x, 10, mu2)
p /= np.sum(p)
b = deepcopy(y) * 0 + 1

r=1
for n in range(0, 1000):
    try:
        c = ud.cconv(b, p)
        b = ud.safedivide(b * y, c)
        d = np.sum(c - y)/np.sum(y)
        e = np.sum(b)/np.sum(y)
        if n>5:
            r=1-0.9*(e-d)
            p=p**r
        print(n, d, e, r)
    except:
        plt.plot(p/np.amax(p))
        plt.plot(b/np.amax(b))
        plt.show()
        exit()

offset = 0
plt.plot(x, b / np.amax(b) - offset * 2)
plt.plot(x, c - offset)
plt.plot(x, y)
plt.show()
