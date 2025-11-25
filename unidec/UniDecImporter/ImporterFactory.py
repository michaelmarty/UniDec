"""To create importer for any of the 4 types,
Simply use: importer_name = ImporterFactory.create_importer(file_path)
To fetch all scans use importer_name.get_all_scans()
"""
import os
import platform
from unidec.UniDecImporter import SingleScanImporter as SSI
from unidec.UniDecImporter.I2MS.I2MS import I2MSImporter
from unidec.UniDecImporter.MZML.mzML import MZMLImporter
from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter

# Note, it is important that these be listed with raw data formats first and processed data formats later.
# Batch.py will attempt the latter formats if use_converted option is on
recognized_types = [".raw", ".mzxml",".mzml", ".mzml.gz", ".gz", '.txt', '.dat', '.csv', '.npz', '.i2ms', '.dmt','.bin']
                     # ".d", ".wiff"]

if platform.system() == "Windows":
    # try:
    #     from Scripts.Importers.Agilent.AgilentImporter import AgilentImporter
    # except Exception as e:
    #     print("Unable to import AgilentImporter:", e)
    #     recognized_types.remove(".d")

    try:
        from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
    except Exception as e:
        print("Unable to import ThermoImporter:", e)
    try:
        from unidec.UniDecImporter.Waters.Waters import WatersDataImporter
    except Exception as e:
        print("Unable to import WatersDataImporter:", e)

    # try:
    #     from Scripts.Importers.Sciex import SciexImporter
    # except Exception as e:
    #     print("Unable to import SciexImporter:", e)

else:
    print("Not importing Thermo, or Waters importers on non-Windows system")
    # Remove .d and .D from recognized types
    # recognized_types.remove(".d")
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
        elif ending==".mzml" or ending==".mzml.gz" or ending==".gz":
            return MZMLImporter(file_path, **kwargs)
        # elif ending==".d":
        #     return AgilentImporter(file_path, **kwargs)
        # elif ending=='.wiff':
        #     return SciexImporter(file_path, **kwargs)
        elif ending == ".txt" or ending == ".dat" or ending == ".csv" or ending == ".npz" or ending == '.bin':
            return SSI.SingleScanImporter(file_path, **kwargs)
        elif ending == ".dmt" or ending == ".i2ms":
            return I2MSImporter(file_path)
        else:
            print("Unsupported file type:", ending, file_path)
            raise IOError("Unsupported file type")
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
    test = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.mzML"
    importer = ImporterFactory.create_importer(test)
    # print(len(dat))



    exit()


#     testfile = "C:\\Python\\UniDec3\\TestSpectra\\test.RAW"
#     testfile = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\LYZ-F319-2-11-22-P1-A1.D"
#     importer = ImporterFactory.create_importer(testfile)