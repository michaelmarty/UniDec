import ctypes
import os
import numpy as np
from pathlib import Path
import sys
import platform

path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

import matplotlib.pyplot as plt
from unidec.IsoDec.match import MatchedPeak, MatchedCollection, IsoDecConfig
from unidec.IsoDec.plots import plot_pks
from unidec.tools import start_at_iso, datachop

current_path = os.path.dirname(os.path.realpath(__file__))

if platform.system() == "Windows":
    dllname = "isodeclib.dll"
elif platform.system() == "Linux":
    dllname = "isodeclib.so"
else:
    dllname = "isodeclib.dylib"

default_dll_path = start_at_iso(dllname, guess=current_path)

if not default_dll_path:
    print("DLL not found anywhere")

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
    # print(config)
    settings = IDSettings()
    settings.phaseres = int(config.phaseres)
    settings.verbose = int(config.verbose)
    settings.peakwindow = int(config.peakwindow)
    settings.peakthresh = float(config.peakthresh)
    settings.minpeaks = int(config.minpeaks)
    settings.css_thresh = float(config.css_thresh)
    settings.matchtol = float(config.matchtol)
    settings.maxshift = int(config.maxshift)
    settings.mzwindow = (config.mzwindow[0], config.mzwindow[1])
    settings.plusoneintwindow = (config.plusoneintwindow[0], config.plusoneintwindow[1])
    settings.knockdown_rounds = int(config.knockdown_rounds)
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
    def __init__(self, dllpath=None):
        if dllpath is None:
            dllpath = default_dll_path

        modelpath = '\\'.join(dllpath.split("\\")[:-1])

        self.modeldir = modelpath

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
            ctypes.c_char_p,  # Model path
        ]
        self.c_lib.predict_charge.restype = ctypes.c_int

        self.c_lib.process_spectrum.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(MPStruct),
            IDSettings,
            ctypes.c_char_p,
        ]

        self.c_lib.DefaultSettings.argtypes = []
        self.c_lib.DefaultSettings.restype = IDSettings

        self.modeldir = modelpath
        # self.modelpath = ctypes.c_char_p(
        #     os.path.join(self.modeldir, "phase_model_8.bin").encode()
        # )
        # Create null pointer for model path
        self.modelpath = None
        self.config = IsoDecConfig()
        # self.determine_model()

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

    def predict_charge(self, centroids, config=None):

        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        # charge = ctypes.c_int(0)
        # self.modelpath = ctypes.c_char_p(None)

        charge = self.c_lib.predict_charge(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(len(cmz)),
            self.modelpath,
            # ctypes.byref(charge),
        )
        return charge

    def process_spectrum(self, centroids, pks=None, config=None, input_type=None):

        if input_type is None:
            input_type = "Peptide"
        type_c = ctypes.c_char_p(input_type.encode('utf-8'))

        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)

        if config is not None:
            self.config = config
        settings = config_to_settings(self.config)

        elems = (MPStruct * n)()
        matchedpeaks = ctypes.cast(elems, ctypes.POINTER(MPStruct))

        nmatched = self.c_lib.process_spectrum(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(n),
            self.modelpath,
            matchedpeaks,
            settings,
            type_c,
        )

        if pks is None:
            pks = MatchedCollection()
        for p in matchedpeaks[:nmatched]:
            if p.z == 0:
                continue
            pk = MatchedPeak(p.z, p.mz, p.avgmass)
            pk.monoiso = p.monoiso
            pk.peakmass = p.peakmass
            pk.avgmass = p.avgmass

            monoisos = np.array(p.monoisos)
            monoisos = monoisos[monoisos > 0]
            pk.monoisos = monoisos

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
            pk.peakint = p.peakint

            # pk.isodist = fast_calc_averagine_isotope_dist(p.monoiso, p.z)
            pk.isodist = np.transpose((isomz, isodist))
            # pk.isodist[:, 1] = p.peakint

            pk.massdist = np.transpose((isomass, isodist))

            pk.startindex = p.startindex
            pk.endindex = p.endindex

            pks.add_peak(pk)
            pks.add_pk_to_masses(pk, 10)
        return pks

    def determine_model(self, default=True):
        if default:
            self.modelpath = None
        else:
            print("Model Directory:", self.modeldir)
            if self.config.phaseres == 4:
                self.modelpath = ctypes.c_char_p(
                    os.path.join(self.modeldir, "phase_model_4.bin").encode('utf-8')
                )
            elif self.config.phaseres == 8:

                self.modelpath = ctypes.c_char_p(
                    os.path.join(self.modeldir, "phase_model_8.bin").encode('utf-8')
                )
            else:
                print("Invalid phase resolution.", self.config.phaseres)
                raise ValueError("Invalid phase resolution.")

            if self.config.verbose:
                print(
                    "Running C code with phaseres of:",
                    self.config.phaseres,
                )


if __name__ == "__main__":
    eng = IsoDecWrapper()

    eng.config.phaseres = 4
    eng.encode(example)
    print(eng.predict_charge(example))
    dat = eng.process_spectrum(example)

    print(dat)
    exit()
    # filepath = "C:\\Data\\IsoNN\\test2.txt"
    filepath = "Z:\\Group Share\\JGP\\js8b05641_si_001\\etd_spectrum.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    spectrum2 = datachop(spectrum, 891.195, 893.757)


    e1 = eng.encode(spectrum2)
    print("Encoded:", e1.shape)

    # z = eng.predict_charge(spectrum2)
    # print("Predicted Charge:", z)

    # exit()
    pks = eng.process_spectrum(spectrum2)
    print(len(pks.peaks))
    # exit()
    plot_pks(pks, centroids=spectrum, show=True, title="C Interface")
    plt.show()
