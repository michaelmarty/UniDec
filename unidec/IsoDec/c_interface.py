import ctypes
import os
import numpy as np
from unidec.IsoDec.match import MatchedPeak, MatchedCollection, IsoDecConfig
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

isodist = ctypes.c_float * 64
matchedinds = ctypes.c_int * 32
#print(isodist)
class MPStruct(ctypes.Structure):
    _fields_ = [('mz', ctypes.c_float),
                ('z', ctypes.c_int),
                ('monoiso', ctypes.c_float),
                ('peakmass', ctypes.c_float),
                ('avgmass', ctypes.c_float),
                ('area', ctypes.c_float),
                ('peakint', ctypes.c_float),
                ('matchedindsiso', ctypes.c_int * 64),
                ('matchedindsexp', ctypes.c_int * 64),
                ('isomz', isodist),
                ('isodist', isodist),
                ('isomass', isodist),
                ('monoisos', ctypes.c_float * 16),
                ('startindex', ctypes.c_int),
                ('endindex', ctypes.c_int),
                ]

class IDSettings(ctypes.Structure):
    _fields_ = [('verbose', ctypes.c_int),
                ('peakwindow', ctypes.c_int),
                ('peakthresh', ctypes.c_float),
                ('minpeaks', ctypes.c_int),
                ('minmatchper', ctypes.c_float),
                ('css_thresh', ctypes.c_float),
                ('matchtol', ctypes.c_float),
                ('maxshift', ctypes.c_int),
                ('mzwindow', ctypes.c_float * 2),
                ('plusoneintwindow', ctypes.c_float * 2),
                ('knockdown_rounds', ctypes.c_int),
                ('min_score_diff', ctypes.c_float),
                ('noisefilter', ctypes.c_float),
                ('minareacovered', ctypes.c_float),
                ('isolength', ctypes.c_int),
                ('mass_diff_c', ctypes.c_float),
                ('adductmass', ctypes.c_float),
                ('minusoneaszero', ctypes.c_int),
                ]


def config_to_settings(config):
    settings = IDSettings()
    settings.verbose = config.verbose
    settings.peakwindow = config.peakwindow
    settings.peakthresh = config.peakthresh
    settings.minpeaks = config.minpeaks
    settings.minmatchper = config.minmatchper
    settings.css_thresh = config.css_thresh
    settings.matchtol = config.matchtol
    settings.maxshift = config.maxshift
    settings.mzwindow = (config.mzwindow[0], config.mzwindow[1])
    settings.plusoneintwindow = (config.plusoneintwindow[0], config.plusoneintwindow[1])
    settings.knockdown_rounds = config.knockdown_rounds
    settings.min_score_diff = config.min_score_diff
    settings.noisefilter = config.noisefilter
    settings.minareacovered = config.minareacovered
    settings.isolength = 64
    settings.mass_diff_c = config.mass_diff_c
    settings.adductmass = config.adductmass
    settings.minusoneaszero = 1
    return settings


class IsoDecWrapper:
    def __init__(self, dllpath=dllpath, modelpath=b"C:\\Python\\UniDec3\\unidec\\IsoDec\\phase_model_2.bin"):
        # void process_spectrum(const double* cmz, const double* cint, const int n, const char* fname, int* charge)
        self.c_lib = ctypes.CDLL(dllpath)
        self.c_lib.predict_charge.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
                                              ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]

        self.c_lib.process_spectrum.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
                                                ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(MPStruct)]
        self.modelpath = ctypes.c_char_p(modelpath)
        self.config = IsoDecConfig()

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

    def process_spectrum(self, centroids, pks=None, config=None):
        #print("Running C Interface")
        cmz = centroids[:, 0].astype(np.float32)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)
        # print(n)

        if config is not None:
            settings = config_to_settings(config)
        else:
            config = self.config
            settings = config_to_settings(config)

        elems = (MPStruct * n)()
        matchedpeaks = ctypes.cast(elems, ctypes.POINTER(MPStruct))
        nmatched = self.c_lib.process_spectrum(cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                               cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n, self.modelpath,
                                               matchedpeaks, settings)
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

            monoisos = np.array(p.monoisos)
            monoisos = monoisos[monoisos > 0]
            pk.monoisos = monoisos
            #print(monoisos)

            if config is not None:
                pk.scan = config.activescan
                pk.ms_order = config.activescanorder
                pk.rt = config.activescanrt

            isodist = np.array(p.isodist)
            isomz = np.array(p.isomz)
            isomass = np.array(p.isomass)

            b1 = isodist > np.amax(isodist) * 0.001
            isodist = isodist[b1]
            isomz = isomz[b1]
            isomass = isomass[b1]
            pk.matchedintensity = np.sum(isodist)
            pk.isodist = np.transpose((isomz, isodist))
            pk.massdist = np.transpose((isomass, isodist))

            pk.startindex = p.startindex
            pk.endindex = p.endindex

            pk.isodist = fast_calc_averagine_isotope_dist(p.monoiso, p.z)
            #pk.isodist[:,0] /= float(p.z)
            pk.isodist[:,1] *= p.peakint
            pks.add_peak(pk)
            pks.add_pk_to_masses(pk, 10)
        return pks


if __name__ == "__main__":
    filepath = "C:\\Data\\IsoNN\\test2.txt"
    #filepath = "C:\\Data\\IsoNN\\test9.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    #spectrum = ud.datachop(spectrum, 702, 708)

    wrapper = IsoDecWrapper()
    #wrapper.process_spectrum(spectrum)
    pks = wrapper.process_spectrum(spectrum)
    print(len(pks.peaks))
    exit()
    plot_pks(pks, centroids=spectrum, show=True, title="C Interface")
