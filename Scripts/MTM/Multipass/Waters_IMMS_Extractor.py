import numpy as np
import matplotlib.pyplot as plt

from unidec.modules.waters_importer import WatersImporter

if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw"

    wi = WatersImporter.WatersDataImporter(path)

    data = wi.readerMS.CombineDriftScan(wi.function, [1,2,3,4])

    print(len(data))
