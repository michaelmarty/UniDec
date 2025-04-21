from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.UniDecImporter.MZML.mzML import *

os.chdir("D:\\Johnny\\IsoDec\\TestData")

if False:
    file2 = "D:\\Johnny\\GitClones\\mzLib\\mzLib\\Test\\DataFiles\\GUACUG_NegativeMode_Sliced.mzML"
    reader = ImporterFactory.create_importer(file2)
    spectra = []
    for i, spectrum in enumerate(reader.msrun):
        spectra.append(spectrum)

    spectrum = get_data_from_spectrum(spectra[0])
    os.chdir("D:\\Johnny\\IsoDec\\TestData")
    np.savetxt("GUACUG_NegativeMode_Sliced.txt", spectrum)
else:
    spectrum = np.loadtxt("GUACUG_NegativeMode_Sliced.txt")

mzrange = [465,470]
mzrange = None

if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

eng = IsoDecEngine(phaseres=8, use_wrapper=False)
eng.config.knockdown_rounds = 5
eng.config.verbose = True
eng.config.adductmass = -1.007276467
eng.config.matchtol = 20
eng.batch_process_spectrum(spectrum, centroided=True)

print("Total peaks: ", len(eng.pks.peaks))


plot_pks(eng.pks, centroids=spectrum, ccolor="k", forcecolor="b", show=False, plotmass=False, nocentroids=False)
plt.show()