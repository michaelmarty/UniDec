import numpy as np
from unidec import modules as ft
import os
import matplotlib.pyplot as plt

fname="extracts.txt"
dname="C:\\Users\Lab.DESKTOP-HJRL65I\Desktop\Deseree's Timing Tests"

data=np.loadtxt(os.path.join(dname,fname))

cv=np.arange(0,100,5)

fits, fitdat=ft.sig_fit(cv,data[:,0])

print(fits[0])

plt.plot(cv,data[:,0])
plt.plot(cv,fitdat)
plt.show()