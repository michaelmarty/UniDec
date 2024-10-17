import tools as ud
import os
import unidec.modules.data_reader as dr
from unidec.UniDecImporter.Importer import ImporterFactory

dr.register()

agilent_files = ["C:\\Data\\DataExamples\\BSA 1mgmL 0.5uL NativeFIA.d",
                 # "C:\\Data\\DataExamples\\AgilentData\\LYZ-F319-2-11-22-P1-A1.D",
                 "C:\\Data\\DataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d",
                 "C:\\Data\\DataExamples\\BSA 10ug FIA 100mM NaOAc flowGrad.d"]

agilent_validation = [11782, 45267, 21124]

for i, f in enumerate(agilent_files):
    print(f)
    file = os.path.abspath(f)
    importer = ImporterFactory.create_importer(file)
    data = importer.get_data()
    l = len(data)
    print("Data Length:", l)
    v = agilent_validation[i]
    if v != l:
        raise Exception("Data length does not match expected value.")

