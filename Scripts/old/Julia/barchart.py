import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib

n_groups = 6
Control_mean = (100, 77.65, 41.287, 19.3267, 10.47, 4.313)
controlStd = (0, 3.189, 4.969, 5.388, 4.060, 2.111)
Imidazole_means = (100, 67.95, 30.77, 11.99, 4.06, 1.15)
imidazolestd = (0, 16.62, 8.358, 1.914, 1.526, 0.0829)
Isopropylimidazole_means = (100, 63.8, 31.47, 7.656, 1.99999, 3.31)
isopropylstd = (0, 7.422, 4.439, 3.15, 0.839, 1.22)
#TMAO_mean = (100, 72.35, 47.03, 30.43, 17.713, 12.15)
#TMAOstd = (0, 12.66, 9.97, 6.299, 4.97, 3.82)

path = "Z:\\Group Share\\Julia\\Imidazole analogues\\190818\\"
os.chdir(path)

yerr = np.loadtxt("stddev.txt")

colors = ["k", "#808080", "b", "r", "c", "m", "y", "g"]
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12



ind = np.arange(n_groups)


space = 0.1


width = (1.5-space)/(n_groups)
plt.bar(ind+width, Control_mean, width, label='Control', color="darkslategrey", yerr=controlStd)
plt.bar(ind+width*2, Imidazole_means, width, label= 'Imidazole', color="aqua", yerr=imidazolestd)
plt.bar(ind+width*3, Isopropylimidazole_means, width, label='2-Isopropylimidazole', color="darkcyan", yerr=isopropylstd)
#plt.bar(ind+width*4, TMAO_mean, width, label='TMAO', color='turquoise', yerr=TMAOstd)
plt.ylabel('Relative Intensity')
plt.xlabel('# of Lipids Bound')
plt.xticks(ind + width / 2, (0, 1, 2, 3, 4, 5))
plt.legend(loc='best')
plt.show()