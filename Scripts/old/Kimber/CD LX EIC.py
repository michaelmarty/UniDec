import os
import glob
import unidec.tools as ud
from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 24
rcParams['font.weight'] = 700
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 700
rcParams['axes.titlesize'] = 16


topdir = "Z:\\Group Share\\Dash\\Data\\Synapt Data\\LXMS\\Test"
files = glob.glob(f'{topdir}/*.raw', recursive=False)

output = []
print(files)

mass1 = 678.384
mass2 = 760.447


mtol = 0.05

outfile = "eic.npz"

for path in files[:]:
    print(path)
    os.chdir(path)

    if not os.path.isfile(outfile):
        try:
             d = WDI(path)
        except Exception as e:
            print("File Failed", path)
            print(e)
            continue
        tic1 = d.get_eic(mass1, mtol)
        tic2 = d.get_eic(mass2, mtol)
        np.savez_compressed(outfile, array1=tic1, array2=tic2)

    path2 = os.path.join(path, outfile)
    print(path2)

    data = np.load(path2)
    tic_1 = data['array1']
    tic_2 = data['array2']
    area1 = ud.integrate(tic_1, 5.25, 6.25)[0]
    area2 = ud.integrate(tic_2, 7.5, 8.5)[0]
    relarea1 = (area1 / (area1 + area2)) * 100
    relarea2 = (area2/ (area1 + area2)) * 100
    print(relarea1, relarea2, area1, area2)
    output.append([path, str(relarea1), str(relarea2), str(area1), str(area2)])

    plt.plot(tic_1[:, 0], tic_1[:, 1], label=str(path) + "DMPC", color = 'green')
    plt.plot(tic_2[:, 0], tic_2[:, 1], label="760.447")


os.chdir(topdir)
np.savetxt("AutoCDLXOG.csv", output, delimiter=",", fmt='%s')
plt.xlim(5, 10)
plt.legend()
plt.show()

