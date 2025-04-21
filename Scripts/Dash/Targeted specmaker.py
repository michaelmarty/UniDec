import unidec.tools as ud
from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
import numpy as np
import matplotlib

def correct_mass(data, exact, window=1):
    d2 = ud.datachop(data, exact - window, exact + window)
    maxval = d2[np.argmax(d2[:, 1]), 0]
    data[:, 0] = data[:, 0] + exact - maxval
    return data


matplotlib.use('WxAgg')

path="Z:\Group Share\Dash\Daptomycin\Control\DJR220501 POPC 18to1 Con3.raw"
monomass = 1619.710
mz = (monomass + 2 * ud.protmass) / 2
mzox = mz + ud.oxmass / 2
window = 0.08
output = []
print(output)


outfile = "spectrum2.npz"


d = WDI(path)
data = d.get_data(time_range=[6.5, 8.5])
np.savez_compressed(outfile, data=data)

#PeakIntegration
forwint=0.08
backint=0.08
mzint = ud.integrate(data, mz - backint, mz + forwint)[0]
oxint = ud.integrate(data, mzox - backint, mzox + forwint)[0]
relox = (oxint / mzint)
print(relox, mzint, oxint)




