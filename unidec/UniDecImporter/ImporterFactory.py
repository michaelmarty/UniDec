"""To create importer for any of the 4 types,
Simply use: importer_name = ImporterFactory.create_importer(file_path)
"""
import os
import platform
import numpy as np

from unidec.UniDecImporter import SingleScanImporter as SSI
from unidec.UniDecImporter.MZML.mzML import MZMLImporter
from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter

# Note, it is important that these be listed with raw data formats first and processed data formats later.
# Batch.py will attempt the latter formats if use_converted option is on
recognized_types = [".raw", ".d", ".mzxml",".mzml", ".mzml.gz", ".gz", '.txt', '.dat', '.csv', '.npz']

if platform.system() == "Windows":
    try:
        from unidec.UniDecImporter.Agilent.AgilentImporter import AgilentImporter
    except Exception as e:
        print("Unable to import AgilentImporter:", e)
        # Remove .d and .D from recognized types
        recognized_types.remove(".d")

    try:
        from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
    except Exception as e:
        print("Unable to import ThermoImporter:", e)
    try:
        from unidec.UniDecImporter.Waters.Waters import WatersDataImporter
    except Exception as e:
        print("Unable to import WatersDataImporter:", e)
else:
    print("Not importing Agilent, Thermo, or Waters importers on non-Windows system")
    # Remove .d and .D from recognized types
    recognized_types.remove(".d")
    # Remove Raw from recognized types
    recognized_types.remove(".raw")


class ImporterFactory:
    def __init__(self):
         self.recognized_file_types = recognized_types

    @staticmethod
    def create_importer(file_path, **kwargs):
        ending = os.path.splitext(file_path)[1]
        ending = ending.lower()
        if ending == ".raw":
            if os.path.isdir(file_path):
                return WatersDataImporter(file_path, **kwargs)
            return ThermoImporter(file_path, **kwargs)
        elif ending==".mzxml":
            return MZXMLImporter(file_path, **kwargs)
        elif ending==".mzml":
            return MZMLImporter(file_path, **kwargs)
        elif ending==".d":
            return AgilentImporter(file_path, **kwargs)
        #Things to introduce in future: .Wiff (Sciex Data)
        elif ending == ".txt" or ending == ".dat" or ending == ".csv" or ending == ".npz":
            return SSI.SingleScanImporter(file_path, **kwargs)
        else:
            print("Unsupported file type:", ending, file_path)
            return None



def get_polarity(path):
    try:
        importer = ImporterFactory.create_importer(path)
        polarity = importer.get_polarity()
    except:
        polarity = "Positive"
    return polarity


if __name__ == "__main__":
    test = u"C:\\Python\\UniDec3\\TestSpectra\\test.raw"
    importer = ImporterFactory.create_importer(test)
    dat = importer.get_data()
    print(len(dat))

    dat = importer.get_data()
    print(len(dat))

    import matplotlib.pyplot as plt

    plt.plot(dat[:, 0], dat[:, 1])
    plt.show()

    exit()


#     testfile = "C:\\Python\\UniDec3\\TestSpectra\\test.RAW"
#     testfile = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\LYZ-F319-2-11-22-P1-A1.D"
#     importer = ImporterFactory.create_importer(testfile)