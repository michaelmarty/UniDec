import sqlite3
import numpy as np


class I2MSImporter:
    def __init__(self, file):
        conn = sqlite3.connect(file)
        cursor = conn.cursor()
        cursor.execute('SELECT * from main.Ion')
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

    def grab_data(self):
        slopes = self.data[:, self.slopekey]
        mz = self.data[:, self.mzkey]
        return np.transpose([mz, slopes])


if __name__ == "__main__":
    file = "Z:\\Group Share\\Marius Kostelic\\STORI\\20220123_MK_BSA_STORI_5_data_2022-02-03-09-10-47.i2MS"
    file = "Z:\\Group Share\\Marius Kostelic\\Baker Lab AAVs\\Replicates STORI for Jack\\AAV2_E3_F3_STORI.dmt"

    i2ms = I2MSImporter(file)
    print(i2ms.keys)
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
