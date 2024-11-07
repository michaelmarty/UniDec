"""To create importer for any of the 4 types,
Simply use: importer_name = ImporterFactory.create_importer(file_path)
"""
import os

recognized_types = [".raw", ".RAW", ".mzxml", ".MZXML", ".mzml", ".MZML", ".d", ".D", ".gz"]

class ImporterFactory:


    def __init__(self):
         self.recognized_file_types = recognized_types

    @staticmethod
    def create_importer(file_path, **kwargs):
        if file_path.endswith(".raw") or file_path.endswith(".RAW"):
            if os.path.isdir(file_path):
                from unidec.UniDecImporter.Waters.Waters import WatersDataImporter
                return WatersDataImporter(file_path, **kwargs)
            from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
            return ThermoImporter(file_path, **kwargs)
        elif file_path.endswith(".mzXML") or file_path.endswith("MZXML") or file_path.endswith("mzxml"):
            from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter
            return MZXMLImporter(file_path, **kwargs)
        elif (file_path.endswith(".mzML") or file_path.endswith(".gz") or file_path.endswith(".MZML") or
              file_path.endswith(".mzml") or file_path.endswith("GZ") or file_path.endswith("mzML")):
            from unidec.UniDecImporter.MZML.mzML import MZMLImporter
            MZMLI = MZMLImporter(file_path, **kwargs)
            return MZMLI
            #return MZMLImporter(file_path, **kwargs)
        elif file_path.endswith(".d") or file_path.endswith(".D"):
            from unidec.UniDecImporter.Agilent.AgilentImporter import AgilentImporter
            return AgilentImporter(file_path, **kwargs)
        #Things to introduce in future: .Wiff (Sciex Data)
        else:
            print("Unsupported file type")
            return None
