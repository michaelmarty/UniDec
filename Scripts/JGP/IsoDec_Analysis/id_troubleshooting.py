from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.match import *
from Scripts.MTM.proteomicstools import *
import os
import time
from UniDecImporter.ImporterFactory import ImporterFactory

file = "Z:\\Group Share\\JGP\\PXD019247\\20170105_L_MaD_ColC4_Ecoli20161108-10_01.raw"

scans = [979]

eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

eng.config.knockdown_rounds = 10

eng.process_file(file, scans)
reader = ImporterFactory.create_importer(file)
spectrum = reader.grab_centroid_data(scans[0])


plot_pks(eng.pks, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=False)
plt.show()