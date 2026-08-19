import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.tools import get_autocorr_ratio, is_centroided
from UniDecImporter.ImporterFactory import ImporterFactory
import scipy.signal as signal

# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

reader = ImporterFactory.create_importer(file2)

centroided = reader.grab_centroid_data(1)
spectrum = reader.get_single_scan(1)

print(len(centroided))
print(len(spectrum))





print(get_autocorr_ratio(centroided))
print(get_autocorr_ratio(spectrum))
