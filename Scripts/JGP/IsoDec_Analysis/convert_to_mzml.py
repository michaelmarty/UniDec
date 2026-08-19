import os
from UniDecImporter.ImporterFactory import ImporterFactory
import numpy as np
import pyopenms
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control\\23_04_21_PEPPI_1B_A.raw"
os.chdir(os.path.dirname(file2))
reader = ImporterFactory.create_importer(file2)
spectrum = reader.get_single_scan(3382)
mz = spectrum[:,0]
intensities = spectrum[:,1]

intensities_zerosremoved = intensities[intensities > 0]
log_intensities_zerosremoved = np.log10(intensities_zerosremoved)

plt.hist(log_intensities_zerosremoved, bins=250)
plt.show()
exit()

plt.plot(mz, intensities)
plt.show()

exp = pyopenms.MSExperiment()

spec = pyopenms.MSSpectrum()
spec.set_peaks((mz,intensities))
spec.set_type(pyopenms.MSSpectrum.Type_PEAKS)
exp.addSpectrum(spec)

pyopenms.MzMLFile().store("uvpd_spectrum.mzML", exp)