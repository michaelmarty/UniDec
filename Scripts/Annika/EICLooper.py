import numpy as np
import pandas as pd
import unidec.tools as ud
import unidec.modules.waters_importer.WatersImporter as WDI
import os


dir = "Z:\\Group Share\\Annika\\APCI Data\\MSMS_ToFMRM Tests\\MS_MSMS_ToFMRM"
os.chdir(dir)

paths = ud.match_dirs_recursive(dir, ".raw")

cmass = 369.2709
startime = 5.4
endtime = 5.8

ints = []
type = []
files = []
concs = []

for i, p in enumerate(paths):
    p = os.path.split(p)[1]
    if "uM" in p:
        # find index of uM
        index = p.find("uM")
        index2 = p[:index].rfind("_")
        conc = p[index2+1:index]
        l = p[:-3] + "_"
    else:
        continue

    if "ToFMRM" in p:
        importer = WDI.WatersDataImporter(p)
        data = importer.get_eic(cmass, tolerance=0.05)
        i1 = ud.integrate(data, startime, endtime)[0]
        ints.append(i1)
        print(i, p, i1)
        files.append(p)
        type.append("ToFMRM")
        concs.append(conc)
        np.savetxt(l + "ToFMRM.txt", data)
    else:
        importer = WDI.WatersDataImporter(p)
        data = importer.get_eic(cmass, tolerance=0.05)
        i1 = ud.integrate(data, startime, endtime)[0]
        ints.append(i1)
        files.append(p)
        type.append("MS")
        concs.append(conc)
        np.savetxt(l + "MS.txt", data)

        importer2 = WDI.WatersDataImporter(p, function=1)
        data2 = importer2.get_eic(cmass, tolerance=0.05)
        i2 = ud.integrate(data2, startime, endtime)[0]
        ints.append(i2)
        files.append(p)
        type.append("MSMS")
        concs.append(conc)
        np.savetxt(l + "MSMS.txt", data2)
        print(i, p, i1, i2)

# create empty data frame
df = pd.DataFrame()
df["Path"] = files
df["Type"] = type
df["Concentration"] = concs
df["Intensity"] = ints
df.to_excel("EIC_Intensities.xlsx", index=False)




