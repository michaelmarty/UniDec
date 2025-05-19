from unidec.IsoDec.plots import plot_pks
from unidec.UniDecImporter import ImporterFactory
from unidec.IsoDec.runtime import IsoDecRuntime
import matplotlib.pyplot as plt
import unidec.tools as ud

filepath = r"Z:\Group Share\JGP\JPST001885\UVPD_Test\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"

importer = ImporterFactory.create_importer(filepath)
spectrum = importer.grab_centroid_data(1)

mz_range = [1000, 1020]
spectrum = ud.datachop(spectrum, mz_range[0], mz_range[1])

eng = IsoDecRuntime(phaseres=8)
eng.config.knockdown_rounds = 5

eng.batch_process_spectrum(spectrum, centroided=True)
print("N Peaks:", len(eng.pks.peaks))

plot_pks(eng.pks, centroids=spectrum, ccolor='k', forcecolor='b', plotmass=True)
plt.show()