import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud

dir = "Z:\\Group Share\\James Rohrbough\\FPOP\\UBIQUITIN TEST\\20190410\\20190410 JGR DESALT10OX3_0_ms_unidecfiles"

file = "20190410 JGR DESALT10OX3_0_ms_mass.txt"

path = os.path.join(dir, file)

data = np.loadtxt(path)

# plt.plot(data[:, 0], data[:, 1])
# plt.show()

ubi = 8565.
ox = np.arange(0, 4) * 16
masses = ubi + ox

extracts=[]
for i in range(0,4):
    ext = ud.data_extract(data, masses[i], 0)
    extracts.append(ext)

extracts = extracts / np.amax(extracts)

plt.bar(range(0,4), extracts)
plt.xticks(range(0,4),ox)
plt.show()