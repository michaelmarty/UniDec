import os
import glob
import unidec.tools as ud
from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
import matplotlib.pyplot as plt
import numpy as np

topdir = "Z:\\Group Share\\Dash\\Data\\Synapt Data\\PLIMB Data February"
files = glob.glob(f'{topdir}/*.raw', recursive=False)

print(files)

mass1 = 817.97
mass2 = 816.95
mass3 = 310.9

masslist = [mass1, mass2, mass3]

mtol = 0.1

for path in files[:]:
    print(path)
    os.chdir(path)
    try:
        d = WDI(path)
    except Exception as e:
        print("File Failed", path)
        print(e)
        continue



    tic = d.get_eic(mass, mtol)
    area1 = ud.integrate(tic, mass1)
    area2 = ud.integrate(tic, mass2)
    area3 = ud.integrate(tic, mass3)
    relarea1 = (area1/(area3+area1))*100
    relarea2 = (area2/(area3+area2))*100
    print(relarea1, relarea2)

    np.savetxt("eic.txt", tic)

    plt.plot(tic[:,0], tic[:,1])

plt.show()
