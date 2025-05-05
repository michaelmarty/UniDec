from unidec.IsoDec.plots import plot_pks
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.UniDecImporter import ImporterFactory
import unidec.tools as ud
import matplotlib.pyplot as plt

file = r"Z:\Group Share\JGP\js8b05641_si_001\1500_scans_200K_16 fills-qb1.raw" ##replace this with your raw file
importer = ImporterFactory.create_importer(file)
spectrum = importer.grab_centroid_data(1) ##Ensure that you're loading the scan you want.

##This is optional if you only want to look at a particular m/z range in the spectrum
mz_range = [700,750]
spectrum = ud.datachop(spectrum, mz_range[0], mz_range[1])

eng = IsoDecRuntime(phaseres=8)
eng.batch_process_spectrum(spectrum)

plot_pks(eng.pks, centroids=spectrum, ccolor='k', forcecolor='g',  plotmass=False)
plt.show()