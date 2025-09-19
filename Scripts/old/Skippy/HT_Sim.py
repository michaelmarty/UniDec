import numpy as np
import unidec.modules.HTEng as HTEng
import matplotlib.pyplot as plt
import unidec.tools as ud

eng = HTEng.HTEng()

# Simulate a gaussian distribution
xdata = np.arange(0, 20, 0.01)
data = ud.ndis(xdata, 3, 0.25)
data = np.transpose(np.vstack((xdata, data)))

# Engine
eng = HTEng.HTEng()
eng.fulltime = xdata
eng.fullscans = np.arange(len(xdata))

eng.config.htbit = -5
eng.setup_ht()

simkernel = eng.htkernel
b1 = simkernel < 0
simkernel[b1] = 0

simdata = np.abs(np.fft.ifft(np.fft.fft(data[:, 1]) * np.fft.fft(simkernel)))
simdata += np.random.normal(0, 0.0, len(simdata))

decon, dat = eng.htdecon(simdata)

plt.plot(data[:, 0], data[:, 1], color="k")
plt.plot(eng.fulltime, simdata, color="y")
plt.plot(eng.fulltime, decon, color="r")
plt.show()