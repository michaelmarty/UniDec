import numpy as np
import matplotlib.pyplot as plt


x = [0.1, 0.2, 0.4, 0.8, 1.6, 6.4]#response curve lipids
y = [0.1823324, 0.307781608, 0.634115738, 0.74174, 1.696481412, 5.950454642]#response curve lipids

#x = [0.270157444, 0.58492946, 0.818901244, 1.637802489, 3.275604977, 6.551209955, 13.10241991]#response curve NDs
#y = [0.169294126, 0.588501864, 0.702237948, 1.410234373, 2.511525356, 5.512220105, 9.184237261]#response curve NDs

coef = np.polyfit(x, y, 1)
poly1d_fn = np.poly1d(coef)
print(coef)

plt.plot(x, y, 'ro', x, poly1d_fn(x), '--k')
plt.plot(x, x, color='b', linestyle="--")
plt.xlim(0, 7)
plt.ylim(0, 6)
plt.xlabel('Prediction', size=20)
plt.ylabel('Actual', size=20)
plt.title("Response Curve", size=20)

plt.show()

#fig.savefig("response lipids.png")