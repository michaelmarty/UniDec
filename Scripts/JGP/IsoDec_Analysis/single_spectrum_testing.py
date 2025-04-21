from Scripts.MTM.proteomicstools import *
from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pyperclip
import pandas as pd
from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.UniDecImporter.MZML.mzML import *
from pyteomics import mgf

file = "Z:\\Group Share\\JGP\\PXD019247\\20170105_L_MaD_ColC4_Ecoli20161116-10_01.raw"

scan = 1832

os.chdir(os.path.dirname(file))

reader = ImporterFactory.create_importer(file)

spectrum = reader.grab_centroid_data(scan)

eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

eng.config.knockdown_rounds = 5
eng.config.css_thresh = 0.70
eng.predmode = 0

eng.batch_process_spectrum(spectrum, centroided=True)

eng.pks.merge_missed_monoisotopics(ppm_tolerance=20, max_mm=3)
for p in eng.pks.peaks:
    print("Peak: ",p.mz, p.monoisos, p.z)

seq = "AKGIREKIK[42.010565]LVSSAGTGHFYTTTKNKRTKPEKLELKKFDPVVRQHVIYKEAKIK"
searchTD(eng.pks.peaks, seq, ["b", "y", "c","z"], "ppm", 20.0, verbose=True)

#Get the peaks with a matched ion as a collection
matched = [p for p in eng.pks.peaks if len(p.matchedions) > 0]
matched_mc = MatchedCollection()
matched_mc.add_peaks(matched)

unmatched = [p for p in eng.pks.peaks if len(p.matchedions) == 0]
unmatched_mc = MatchedCollection()
unmatched_mc.add_peaks(unmatched)


plot_pks(matched_mc, centroids=spectrum, ccolor="k", forcecolor="b", show=True, plotmass=False, nocentroids=True, labelz=True, labelpeaks=True)
plot_pks(unmatched_mc, centroids=spectrum, ccolor="k", forcecolor="r", show=True, plotmass=False, nocentroids=True, labelz=True, labelpeaks=True)
plt.show()