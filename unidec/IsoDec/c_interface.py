import ctypes
import os
import numpy as np
from unidec.IsoDec.match import MatchedPeak, MatchedCollection, IsoDecConfig
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist
from unidec.IsoDec.plots import *
import unidec.tools as ud
#from unidec.IsoDec.encoding import encode_phase

example = np.array(
    [
        [5.66785531e02, 1.47770838e06],
        [5.67057354e02, 1.54980838e06],
        [5.67507468e02, 5.21600520e07],
        [5.67708173e02, 8.35557760e07],
        [5.67908401e02, 7.28264240e07],
        [5.68060254e02, 1.87337225e06],
        [5.68108674e02, 4.35435520e07],
        [5.68239256e02, 3.88155375e06],
        [5.68309390e02, 2.05468060e07],
        [5.68509951e02, 7.18109250e06],
        [5.68707871e02, 2.30373500e06],
        [5.69150563e02, 1.57598062e06],
        [5.69243121e02, 1.96390440e07],
        [5.69334393e02, 6.82677120e07],
        [5.69425337e02, 1.22867432e08],
        [5.69516492e02, 1.45702336e08],
        [5.69607541e02, 1.20801936e08],
        [5.69698595e02, 1.06786072e08],
        [5.69789906e02, 6.56232960e07],
        [5.69881208e02, 3.41013880e07],
        [5.69972168e02, 1.70930360e07],
        [5.70063432e02, 9.17621100e06],
        [5.70699369e02, 1.96462650e06],
    ]
)

dllpath = (
    "C:\\Python\\UniDec3\\unidec\\IsoDec\\src\\isodec\\x64\\Release\\isodeclib.dll"
)

isodist = ctypes.c_float * 64
matchedinds = ctypes.c_int * 32


# print(isodist)
class MPStruct(ctypes.Structure):
    _fields_ = [
        ("mz", ctypes.c_float),
        ("z", ctypes.c_int),
        ("monoiso", ctypes.c_float),
        ("peakmass", ctypes.c_float),
        ("avgmass", ctypes.c_float),
        ("area", ctypes.c_float),
        ("peakint", ctypes.c_float),
        ("matchedindsiso", ctypes.c_int * 64),
        ("matchedindsexp", ctypes.c_int * 64),
        ("isomz", isodist),
        ("isodist", isodist),
        ("isomass", isodist),
        ("monoisos", ctypes.c_float * 16),
        ("startindex", ctypes.c_int),
        ("endindex", ctypes.c_int),
        ("score", ctypes.c_float),
        ("realisolength", ctypes.c_int),
    ]


class IDSettings(ctypes.Structure):
    _fields_ = [
        ("phaseres", ctypes.c_int),
        ("verbose", ctypes.c_int),
        ("peakwindow", ctypes.c_int),
        ("peakthresh", ctypes.c_float),
        ("minpeaks", ctypes.c_int),
        ("css_thresh", ctypes.c_float),
        ("matchtol", ctypes.c_float),
        ("maxshift", ctypes.c_int),
        ("mzwindow", ctypes.c_float * 2),
        ("plusoneintwindow", ctypes.c_float * 2),
        ("knockdown_rounds", ctypes.c_int),
        ("min_score_diff", ctypes.c_float),
        ("minareacovered", ctypes.c_float),
        ("isolength", ctypes.c_int),
        ("mass_diff_c", ctypes.c_double),
        ("adductmass", ctypes.c_float),
        ("minusoneaszero", ctypes.c_int),
        ("isotopethreshold", ctypes.c_float),
        ("datathreshold", ctypes.c_float),
        ("zscore_threshold", ctypes.c_float),
    ]


class IDConfig(ctypes.Structure):
    _fields_ = [
        ("verbose", ctypes.c_int),
        ("pres", ctypes.c_int),
        ("maxz", ctypes.c_int),
        ("elen", ctypes.c_int),
        ("l1", ctypes.c_int),
        ("l2", ctypes.c_int),
        ("l3", ctypes.c_int),
        ("l4", ctypes.c_int),
        ("dlen", ctypes.c_int),
    ]


def config_to_settings(config):
    settings = IDSettings()
    settings.phaseres = config.phaseres
    settings.verbose = config.verbose
    settings.peakwindow = config.peakwindow
    settings.peakthresh = config.peakthresh
    settings.minpeaks = config.minpeaks
    settings.css_thresh = config.css_thresh
    settings.matchtol = config.matchtol
    settings.maxshift = config.maxshift
    settings.mzwindow = (config.mzwindow[0], config.mzwindow[1])
    settings.plusoneintwindow = (config.plusoneintwindow[0], config.plusoneintwindow[1])
    settings.knockdown_rounds = config.knockdown_rounds
    settings.min_score_diff = config.min_score_diff
    settings.minareacovered = config.minareacovered
    settings.isolength = 64
    settings.mass_diff_c = config.mass_diff_c
    settings.adductmass = config.adductmass
    settings.minusoneaszero = config.minusoneaszero
    settings.isotopethreshold = config.isotopethreshold
    settings.datathreshold = config.datathreshold
    settings.zscore_threshold = config.zscore_threshold
    return settings


class IsoDecWrapper:
    def __init__(self, dllpath=None, modeldir="C:\\Python\\UniDec3\\unidec\\IsoDec\\"):
        if dllpath is None:
            dllpath = os.path.join(modeldir, "isodeclib.dll")

        self.c_lib = ctypes.CDLL(dllpath)

        self.c_lib.encode.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float),
            IDConfig,
            IDSettings,
        ]

        self.c_lib.predict_charge.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
        ]

        self.c_lib.process_spectrum.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(MPStruct),
        ]
        self.modeldir = modeldir
        self.modelpath = ctypes.c_char_p(
            os.path.join(self.modeldir, "phase_model.bin").encode()
        )
        self.config = IsoDecConfig()

    def encode(self, centroids, maxz=50, phaseres=8, config=None):
        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        elen = maxz * phaseres
        emat = np.zeros(elen).astype(np.float32)
        idconf = IDConfig()
        idconf.pres = phaseres
        idconf.maxz = maxz
        idconf.elen = elen

        if config is not None:
            self.config = config
            settings = config_to_settings(config)
        else:
            config = self.config
            settings = config_to_settings(config)

        settings = config_to_settings(self.config)
        self.c_lib.encode(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(len(cmz)),
            emat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            idconf,
            settings,
        )
        # Convert emat to numpy
        emat = np.ctypeslib.as_array(emat)
        return emat

    def predict_charge(self, centroids):
        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)
        charge = ctypes.c_int(0)
        self.c_lib.predict_charge(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n,
            self.modelpath,
            ctypes.byref(charge),
        )
        return charge.value

    def process_spectrum(self, centroids, pks=None, config=None):
        # print("Running C Interface")
        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)
        # print(n)

        if config is not None:
            self.config = config
            settings = config_to_settings(config)
        else:
            config = self.config
            settings = config_to_settings(config)

        if self.config.phaseres == 4:
            self.modelpath = ctypes.c_char_p(
                os.path.join(self.modeldir, "phase_model_2.bin").encode()
            )
        elif self.config.phaseres == 8:
            self.modelpath = ctypes.c_char_p(
                os.path.join(self.modeldir, "phase_model.bin").encode()
            )
        else:
            print("Invalid phase resolution.", self.config.phaseres)
            raise ValueError("Invalid phase resolution.")

        if config.verbose:
            print(
                "Running C code with",
                str(self.modelpath.value),
                "model and phaseres of:",
                self.config.phaseres,
            )
        elems = (MPStruct * n)()
        matchedpeaks = ctypes.cast(elems, ctypes.POINTER(MPStruct))
        nmatched = self.c_lib.process_spectrum(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n,
            self.modelpath,
            matchedpeaks,
            settings,
        )
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
            # print(monoisos)

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
            # pk.isodist[:,0] /= float(p.z)
            pk.isodist[:, 1] *= p.peakint
            pks.add_peak(pk)
            pks.add_pk_to_masses(pk, 10)
        return pks


if __name__ == "__main__":
    filepath = "C:\\Data\\IsoNN\\test2.txt"
    # filepath = "C:\\Data\\IsoNN\\test9.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    spectrum = ud.datachop(spectrum, 702, 708)
    spectrum = spectrum.astype(np.double)

    wrapper = IsoDecWrapper()
    e1 = wrapper.encode(spectrum)
    exit()
    e2 = encode_phase(spectrum).flatten()

    diff = e1 - e2
    #diff[np.abs(diff)<0.0001] = 0
    print(np.allclose(e1, e2))
    print(diff)
    print(np.amax(np.abs(diff)))
    print("Total:", np.sum(np.abs(diff)))
    exit()
    # wrapper.process_spectrum(spectrum)
    pks = wrapper.process_spectrum(spectrum)
    print(len(pks.peaks))
    exit()
    plot_pks(pks, centroids=spectrum, show=True, title="C Interface")
