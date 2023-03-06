import numpy as np
import os
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLWriter import ImzMLWriter
import h5py
import unidec.tools as ud
from unidec.modules.hdf5_tools import replace_dataset


class Imzml_Reader:
    def __init__(self, path):
        print("Opening iMZML filer:", path)
        self.path = path
        self.dir = os.path.dirname(self.path)
        self.header = os.path.splitext(self.path)[0]
        print(self.header)
        self.msrun = ImzMLParser(path)
        self.coords = self.msrun.coordinates
        self.maxc = np.amax(self.coords, axis=0)
        print(self.maxc)
        self.data = []

    def get_data(self):
        self.data = []
        for idx, (x, y, z) in enumerate(self.coords):
            mzs, intensities = self.msrun.getspectrum(idx)
            self.data.append(np.transpose([mzs, intensities]))
        # self.data = np.array(self.data)
        return self.data

    def write_to_hdf5(self, outfile=None):
        if outfile is None:
            outfile = self.header + ".hdf5"

        hdf = h5py.File(outfile, "a")
        try:
            del hdf["ms_dataset"]
        except:
            pass
        msdataset = hdf.require_group("ms_dataset")

        msdataset.attrs["v1name"] = "xpos"
        msdataset.attrs["v2name"] = "ypos"

        config = hdf.require_group("config")
        config.attrs["metamode"] = -1

        if len(self.data)==0:
            self.get_data()
        num = 0
        for idx, (x, y, z) in enumerate(self.coords):
            data = self.data[idx]
            if not ud.isempty(data):
                group = msdataset.require_group(str(idx))
                replace_dataset(group, "raw_data", data=data)
                group.attrs["Position"] = str(int(x))+"," + str(int(y)) + "," + str(int(z))
                group.attrs["xpos"] = int(x)
                group.attrs["ypos"] = int(y)
                group.attrs["zpos"] = int(z)
                num += 1
                pass

        msdataset.attrs["num"] = num
        hdf.close()

def hdf5_to_imzml(path, outpath):
    from unidec.metaunidec.mudeng import MetaUniDec
    eng = MetaUniDec()
    eng.open(path)

    with ImzMLWriter(outpath) as w:
        for i, s in enumerate(eng.data.spectra):
            mzs = s.massdat[:,0]
            intensities = s.massdat[:,1]
            if len(mzs)<2:
                mzs = [eng.config.masslb, eng.config.massub]
                intensities = [0, 0]
            x = s.attrs['xpos']
            y = s.attrs['ypos']
            z = s.attrs['zpos']
            coords = tuple(np.array([int(x), int(y), int(z)]).astype(int))
            print(coords)
            print(i, coords)
            w.addSpectrum(mzs, intensities, coords)

def imzml_to_hdf5(infile, outfile):
    r = Imzml_Reader(infile)
    data = r.get_data()
    r.write_to_hdf5(outfile)


if __name__ == '__main__':
    path = "C:\\Data\\Imaging\\UoB_RepairedRK_image_200x200um_proc\\RepairedRK_image_200x200um_proc.hdf5"
    outpath = "C:\\Data\\Imaging\\UoB_RepairedRK_image_200x200um_proc\\test.imzml"
    hdf5_to_imzml(path, outpath)

    #exit()
    #path = "C:\\Data\\Imaging\\UoB_RepairedRK_image_200x200um_proc\\RepairedRK_image_200x200um_proc.imzML"
    #r = Imzml_Reader(path)
    #data = r.get_data()
    #r.write_to_hdf5()

    # os.chdir("C:\\Data\\Imaging\\UoBham_nano-DESI_imaging_kidney\\imzML\\")
    # h = "RK_image_200x200um"
    # ImzMLParser("RK_image_200x200um.imzml")
