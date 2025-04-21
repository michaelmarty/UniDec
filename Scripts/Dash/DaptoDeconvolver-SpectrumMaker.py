import os
import glob
import unidec.tools as ud
from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from unidec.engine import UniDec as UD


def correct_mass(data, exact, window=1):
    d2 = ud.datachop(data, exact - window, exact + window)
    maxval = d2[np.argmax(d2[:, 1]), 0]
    data[:, 0] = data[:, 0] + exact - maxval
    return data


def get_vars(path):
    lipids = ["DMPC", "DMPG", "Ecoli", "Free", "DPPC", "DLPC", "POPC", "fdap", "fdapb"]
    ratios = ["3to", "9to", "18to"]
    lipid = ""
    ratio = ""
    for l in lipids:
        if l.lower() in path.lower():
            lipid = l
    for r in ratios:
        if r.lower() in path.lower():
            ratio = r[:-2]
    print(lipid, ratio)
    return lipid, ratio


def run_deconvolution(path, configpath=None):
    eng = UD()
    eng.open_file(path)
    if configpath is not None:
        print("Loading Config:", configpath)
        eng.load_config(configpath)
    eng.run_unidec()
    return eng.data.massdat


matplotlib.use('WxAgg')

topdir = "Z:\Group Share\Dash\Data\Synapt Data\FPOP Data September"
files = glob.glob(f'{topdir}/*.raw', recursive=False)
configpath = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\Config Files\\Regular_Daptomycin_FPOP_Settings.dat"


monomass = 1619.710
mz = (monomass + 2 * ud.protmass) / 2
mzox = mz + ud.oxmass / 2
window = 1
window1 = 0.1
window2 = 0.1
output = []
print(files)


outfile = "spectrum.npz"
for i, path in enumerate(files[:]):
    print(path)
    os.chdir(path)

    if not os.path.isfile(outfile):
        try:
            d = WDI(path)
        except Exception as e:
            print("File Failed", path)
            print(e)
            continue

        data = d.get_data(time_range=[9, 12])
        np.savez_compressed(outfile, data=data)


    path2 = os.path.join(path, outfile)
    print(path2)
    data = run_deconvolution(path2, configpath)

    data = correct_mass(data, mz, window=1)
    data = ud.datachop(data, 805, 830)

    mzint = ud.integrate(data, mz - window, mz + window)[0]
    oxint = ud.integrate(data, mzox - window, mzox + window)[0]
    relox = (oxint / mzint)*100
    print(relox, mzint, oxint)
    if mzint < 10000:
        print("ERROR: ", path)
    lipid, ratio = get_vars(path)
    output.append([path, lipid, ratio, str(relox)])

    plt.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]) - i * 0.05)
    plt.vlines(mz - window, -5, 1)
    plt.vlines(mz + window, -5, 1)
    plt.vlines(mzox - window, -5, 1)
    plt.vlines(mzox + window, -5, 1)
    # plt.xlim(805, 830)

os.chdir(topdir)
np.savetxt("OutputDaptoDeconvolver.csv", output, delimiter=",", fmt='%s')
plt.show()
