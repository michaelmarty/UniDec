"""To create importer for any of the 4 types,
Simply use: importer_name = ImporterFactory.create_importer(file_path)
"""
import os

try:
    from unidec.UniDecImporter.Agilent.AgilentImporter import AgilentImporter
except Exception as e:
    print(e)
    print("Unable to import AgilentImporter")

from unidec.UniDecImporter.MZML.mzML import MZMLImporter
from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter
try:
    from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
except Exception as e:
    print(e)
    print("Unable to import ThermoImporter")
try:
    from unidec.UniDecImporter.Waters.Waters import WatersDataImporter
except Exception as e:
    print(e)
    print("Unable to import WatersDataImporter")

#Should just be case insensitive comparisons when generating importers
recognized_types = [".raw", ".RAW", ".mzxml", ".MZXML", ".mzml", ".MZML", ".mzML", ".d", ".D", ".gz", '.mzXML',  '.mzMl']

class ImporterFactory:
    def __init__(self):
         self.recognized_file_types = recognized_types

    @staticmethod
    def create_importer(file_path, **kwargs):
        if file_path.endswith(".raw") or file_path.endswith(".RAW"):
            if os.path.isdir(file_path):
                return WatersDataImporter(file_path, **kwargs)
            return ThermoImporter(file_path, **kwargs)
        elif file_path.endswith(".mzXML") or file_path.endswith("MZXML") or file_path.endswith("mzxml"):
            return MZXMLImporter(file_path, **kwargs)
        elif (file_path.endswith(".mzML") or file_path.endswith(".gz") or file_path.endswith(".MZML") or
              file_path.endswith(".mzml") or file_path.endswith("GZ") or file_path.endswith("mzML") or file_path.endswith(".mzMl")):
            return MZMLImporter(file_path, **kwargs)
        elif file_path.endswith(".d") or file_path.endswith(".D"):
            return AgilentImporter(file_path, **kwargs)
        #Things to introduce in future: .Wiff (Sciex Data)
        else:
            print("Unsupported file type")
            return None

if __name__ == "__main__":
    print("test")