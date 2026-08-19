from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import plot_pks, gen_plot_fig
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import *
import pyperclip
from Scripts.MTM.proteomicstools import *
from UniDecImporter.ImporterFactory import ImporterFactory
from unidec.IsoDec.encoding import *
from unidec.IsoDec.models import *


# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
os.chdir(os.path.dirname(file2))


eng = IsoDecEngine(phaseres=8, use_wrapper=False)

eng.config.verbose = True
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

eng.config.minareacovered = 0.2
eng.config.knockdown_rounds = 10
eng.config.css_thresh = 0.7
eng.predmode = 4

if False:
    reader = ImporterFactory.create_importer(file2)
    spectrum = reader.grab_centroid_data(1)
    np.savetxt("uvpd_spectrum.txt", spectrum)
else:
    spectrum = np.loadtxt("uvpd_spectrum.txt")

mzrange = [316.5,317.75]
#mzrange = None

if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

#Find the most intense peakmz in the spectrum
maxint = np.amax(spectrum[:,1])
maxindex = np.argmax(spectrum[:,1])
maxmz = spectrum[maxindex,0]

emat = encode_phase(spectrum, maxz=50, phaseres=8)

z, predvec = eng.phasemodel.predict_returnvec(emat)

#order the predvec by descending
order = np.argsort(predvec[0])[::-1]



#Now iterate through the predictions in order to see how that looks
for i in order:
    print("Z:", i, "Prediction:", predvec[0][i])
    matches = eng.get_matches(spectrum, i, maxmz, eng.pks)
    plot_pks(eng.pks, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=True, labelz=True)
    plt.show()

    eng.pks = MatchedCollection()