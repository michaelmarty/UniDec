import os
from UniDecImporter.ImporterFactory import *
import matplotlib.pyplot as plt
import numpy as np

pamass = 255.23242

path = r"Z:\Group Share\Annika\Stellar\20250902\20250902_AS_02_EPL_Targeted_DDA.raw"

im = ImporterFactory()
importer = im.create_importer(path)

localmaxlist = importer.get_mz_localmax(pamass, 1000)
print(localmaxlist)

# Calculate the weighted average of the local maxima positions
if len(localmaxlist) > 0:
    weighted_avg_mz = np.average(localmaxlist[:, 0], weights=localmaxlist[:, 1])
    print("Weighted Average m/z:", weighted_avg_mz)


plt.subplot(221)
plt.scatter(localmaxlist[:,0], localmaxlist[:,1])
plt.subplot(222)
plt.plot(localmaxlist[:,0])
plt.subplot(223)
plt.plot(localmaxlist[:,1])
plt.show()


