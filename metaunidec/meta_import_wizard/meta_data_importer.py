import os
from metaunidec.meta_import_wizard import MetaTagTypes as tt
import metaunidec.mudeng as mudeng
from unidec_modules.unidectools import get_importer
from unidec_modules.mzMLimporter import *

def parse_file(file_path, exp_type='Time', collision=None, dir=None):
    file_name = os.path.split(file_path)[1]
    exedir = dir

    out = {tt.FILE_PATH: file_path,
           tt.FILE_NAME: file_name,
           tt.V1: None,
           tt.V2: None,
           tt.TIME_START: None,
           tt.TIME_END: None,
           tt.TYPE: exp_type,
           tt.SCAN_START: None,
           tt.SCAN_END: None}

    return out


def auto_from_wizard(data, filename, mode):
    # print data
    eng = mudeng.MetaUniDec()
    eng.data.new_file(filename)
    for i, d in enumerate(data):
        v1 = d[1]
        v2 = d[2]
        f=d[0]
        try:
            start = float(d[3])
        except:
            start = None
        try:
            stop = float(d[4])
        except:
            stop = None
        path = d[5]

        if start is None or stop is None:
            eng.data.add_file(path=path)
        else:
            print(start, stop)
            importer = get_importer(path)
            if mode == 1:
                data = importer.get_data(scan_range=(start, stop))
            elif mode == 0:
                data = importer.get_data(time_range=(start, stop))
            eng.data.add_data(data, path)
        eng.data.spectra[-1].var1 = v1
        eng.data.spectra[-1].var2 = v2
        eng.data.spectra[-1].name = f
    eng.data.export_hdf5()


