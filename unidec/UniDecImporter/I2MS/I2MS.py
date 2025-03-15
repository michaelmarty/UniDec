import sqlite3
from collections import Counter

import numpy as np



class I2MSImporter:
    def __init__(self, file):
        conn = sqlite3.connect(file)
        cursor = conn.cursor()
        cursor.execute('SELECT * from main.Ion ORDER BY mz ASC')
        data = cursor.fetchall()
        self.data = np.array(data)
        cursor.execute('PRAGMA table_info(Ion)')
        keys = cursor.fetchall()
        self.keys = np.array(keys)
        self.mzkey = np.where(self.keys[:, 1] == "Mz")[0][0]
        self.scankey = np.where(self.keys[:, 1] == "ScanNumber")[0][0]
        self.slopekey = np.where(self.keys[:, 1] == "Slope")[0][0]
        self.scans = self.data[:,self.scankey]
        self.iitkey = np.where(self.keys[:, 1] == "InverseInjectionTimeSeconds")[0][0]
        self.invinjtime = self.data[:, self.iitkey]
        '''
        try:
            self.invinjtime = np.where(self.keys[:, 1] == "InverseInjectionTimeSeconds")[0][0]
        except:
            self.invinjtime = None'''
        self.cdms_support=True
        self.imms_support=False
        self.chrom_support=True

    def get_all_scans(self, threshold=-1):
        slopes = self.data[:, self.slopekey]
        mz = self.data[:, self.mzkey]

        if threshold >= 0:
            mask = slopes > threshold
            slopes = slopes[mask]
            mz = mz[mask]
        return np.transpose([mz, slopes])

    def get_cdms_data(self):
        raw_dat = self.get_all_scans()
        mz = raw_dat[:, 0]
        intensity = raw_dat[:, 1]
        scans = self.data[:, self.scankey]
        if len(scans) != len(mz):
            scans = np.ones(len(mz))
        it = self.invinjtime
        data_array = np.transpose([mz, intensity, scans, it])
        return data_array

    def get_cdms_data_by_scans(self, scan_range):
        data = self.get_cdms_data()
        mask = np.logical_and(data[:, 2] >= scan_range[0], data[:, 2] <= scan_range[1])
        return data[mask]

    def get_single_scan(self, scan=None):
        res = self.data[self.data[:, self.scankey] == scan]
        return res

    def close(self):
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()


    def get_scan_range(self):
        return [1,  int(max(self.scans))]

if __name__ == "__main__":


    file = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\CDMS\\test_dmt_cdms.dmt"
    i2ms = I2MSImporter(file)
    i2ms.get_all_scans()
    # for i in i2ms.scans:
    #     print(i)


    exit()
    import matplotlib.pyplot as plt

    x = i2ms.data[:,7]
    y = i2ms.data[:,10]

    mzbins = 1000
    zbins = 0.1

    mzrange = [np.floor(np.amin(x)), np.amax(x)/5]
    zrange = [np.floor(np.amin(y)), np.amax(y)]
    mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + mzbins / 2, mzbins)
    # Weird fix to make this axis even is necessary for CuPy fft for some reason...
    if len(mzaxis) % 2 == 1:
        mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + 3 * mzbins / 2, mzbins)
    zaxis = np.arange(zrange[0] - zbins / 2., zrange[1] + zbins / 2, zbins)

    harray, xtab, ytab = np.histogram2d(x, y, [mzaxis, zaxis])
    xtab = xtab[1:] - mzbins / 2.
    ytab = ytab[1:] - zbins / 2.

    harray = np.transpose(harray)
    harray /= np.amax(harray)

    plt.imshow(harray, aspect="auto")
    plt.show()
