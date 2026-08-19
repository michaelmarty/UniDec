from unidec.IsoDec.compare_tools import compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pyperclip
from unidec.IsoDec.match import *
from UniDecImporter.Importer import ImporterFactory

# Set backend to Agg
mpl.use("WxAgg")

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

reader = ImporterFactory.create_importer(file2)
spectrum = reader.grab_centroid_data(1)

spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

mzrange = [1196, 1201]
spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])
spectrum.astype(np.float32)
#spectrum = spectrum[:2]
verbose = True

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.config.datathreshold = 0

print("Pred1:", eng.phase_predictor(spectrum))

#eng.batch_process_spectrum(spectrum, centroided=True)
#mc1 = eng.pks

eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.config.datathreshold = 0

print("Pred2:", eng.phase_predictor(spectrum))

exit()

eng.batch_process_spectrum(spectrum, centroided=True)
mc2 = eng.pks

shared, u1, u2 = compare_matchedcollections(
    mc1, mc2, ppmtol=50, objecttocompare="monoiso", maxshift=3
)


plt.figure()
plot_pks(
    u1,
    forcecolor="r",
    centroids=spectrum,
    ccolor="k",
    show=False,
    plotmass=False,
    nocentroids=True,
)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True)
plot_pks(shared, forcecolor="g", show=False, plotmass=False, nocentroids=True)
plt.show()
exit()
