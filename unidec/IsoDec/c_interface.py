import ctypes
import os
import numpy as np
from unidec.IsoDec.match import MatchedPeak, MatchedCollection
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist
from unidec.IsoDec.plots import *
import unidec.tools as ud

example = np.array([[5.66785531e+02, 1.47770838e+06],
                    [5.67057354e+02, 1.54980838e+06],
                    [5.67507468e+02, 5.21600520e+07],
                    [5.67708173e+02, 8.35557760e+07],
                    [5.67908401e+02, 7.28264240e+07],
                    [5.68060254e+02, 1.87337225e+06],
                    [5.68108674e+02, 4.35435520e+07],
                    [5.68239256e+02, 3.88155375e+06],
                    [5.68309390e+02, 2.05468060e+07],
                    [5.68509951e+02, 7.18109250e+06],
                    [5.68707871e+02, 2.30373500e+06],
                    [5.69150563e+02, 1.57598062e+06],
                    [5.69243121e+02, 1.96390440e+07],
                    [5.69334393e+02, 6.82677120e+07],
                    [5.69425337e+02, 1.22867432e+08],
                    [5.69516492e+02, 1.45702336e+08],
                    [5.69607541e+02, 1.20801936e+08],
                    [5.69698595e+02, 1.06786072e+08],
                    [5.69789906e+02, 6.56232960e+07],
                    [5.69881208e+02, 3.41013880e+07],
                    [5.69972168e+02, 1.70930360e+07],
                    [5.70063432e+02, 9.17621100e+06],
                    [5.70699369e+02, 1.96462650e+06]])

dllpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\src\\isodec\\x64\\Release\\isodeclib.dll"

isodist = ctypes.c_float * 128
#print(isodist)
class MPStruct(ctypes.Structure):
    _fields_ = [('mz', ctypes.c_float),
                ('z', ctypes.c_int),
                ('monoiso', ctypes.c_float),
                ('peakmass', ctypes.c_float),
                ('avgmass', ctypes.c_float),
                ('area', ctypes.c_float),
                ('peakint', ctypes.c_float),
                ('isomz', isodist),
                ('isodist', isodist),
                ('isomass', isodist),
                ('startindex', ctypes.c_int),
                ('endindex', ctypes.c_int),
                ]


class IsoDecWrapper:
    def __init__(self, dllpath=dllpath, modelpath=b"C:\\Python\\UniDec3\\unidec\\IsoDec\\phase_model_2.bin"):
        # void process_spectrum(const double* cmz, const double* cint, const int n, const char* fname, int* charge)
        self.c_lib = ctypes.CDLL(dllpath)
        self.c_lib.predict_charge.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
                                              ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]

        self.c_lib.process_spectrum.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
                                                ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(MPStruct)]
        self.modelpath = ctypes.c_char_p(modelpath)

    def predict_charge(self, centroids):
        cmz = centroids[:, 0].astype(np.float32)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)
        charge = ctypes.c_int(0)
        self.c_lib.predict_charge(cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                  cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n, self.modelpath,
                                  ctypes.byref(charge))
        print(charge.value)
        return charge.value

    def process_spectrum(self, centroids, pks=None, scannumber=None):
        #print("Running C Interface")
        cmz = centroids[:, 0].astype(np.float32)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)
        # print(n)
        elems = (MPStruct * n)()
        matchedpeaks = ctypes.cast(elems, ctypes.POINTER(MPStruct))
        nmatched = self.c_lib.process_spectrum(cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                               cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n, self.modelpath,
                                               matchedpeaks)
        # print(nmatched)
        if pks is None:
            pks = MatchedCollection()
        for p in matchedpeaks[:nmatched]:
            if p.z == 0:
                continue
            pk = MatchedPeak(p.z, p.mz)
            pk.monoiso = p.monoiso
            pk.peakmass = p.peakmass
            pk.avgmass = p.avgmass
            pk.scan = scannumber

            isodist = np.array(p.isodist)
            isomz = np.array(p.isomz)
            isomass = np.array(p.isomass)

            b1 = isodist > np.amax(isodist) * 0.001
            isodist = isodist[b1]
            isomz = isomz[b1]
            isomass = isomass[b1]

            pk.isodist = np.transpose((isomz, isodist))
            pk.massdist = np.transpose((isomass, isodist))

            pk.startindex = p.startindex
            pk.endindex = p.endindex

            #pk.isodist = fast_calc_averagine_isotope_dist(p.monoiso, p.z)
            #pk.isodist[:,0] /= float(p.z)
            #pk.isodist[:,1] *= p.peakint
            pks.add_peak(pk)

        return pks


if __name__ == "__main__":
    filepath = "C:\\Data\\IsoNN\\test2.txt"
    #filepath = "C:\\Data\\IsoNN\\test9.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    #spectrum = ud.datachop(spectrum, 702, 708)

    wrapper = IsoDecWrapper()
    #wrapper.process_spectrum(spectrum)
    pks = wrapper.process_spectrum(spectrum)
    #exit()
    plot_pks(pks, centroids=spectrum, show=True, title="C Interface")
