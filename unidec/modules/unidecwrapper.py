import ctypes
import os
import numpy as np
from pathlib import Path
import sys
import platform
import time

from unidec.tools import start_at_iso
from unidec.engine import UniDec

path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
current_path = os.path.dirname(os.path.realpath(__file__))

default_dll_path = os.path.join(os.path.dirname(current_path), "bin")

if platform.system() == "Windows":
    dllname = "unideclib.dll"
elif platform.system() == "Linux":
    dllname = "unideclib.so"
else:
    dllname = "unideclib.dylib"

default_dll_path = start_at_iso(dllname, guess=default_dll_path)

if not default_dll_path:
    print("DLL not found anywhere")
else:
    print(f"Using DLL at {default_dll_path}")

"""
struct Input {
    float *dataMZ;
    float *dataInt;
    float *testmasses;
    int *nztab;
    float *mtab;
    char *barr;
};
"""
class Input(ctypes.Structure):
    _fields_ = [
        ("dataMZ", ctypes.POINTER(ctypes.c_float)),
        ("dataInt", ctypes.POINTER(ctypes.c_float)),
        ("testmasses", ctypes.POINTER(ctypes.c_float)),
        ("nztab", ctypes.POINTER(ctypes.c_int)),
        ("mtab", ctypes.POINTER(ctypes.c_float)),
        ("barr", ctypes.POINTER(ctypes.c_char))
    ]


def py_to_c_input(dataobj, c_config, py_config):
    """
    Convert a Python data object to a C Input structure.
    :param dataobj: Dictionary containing input data.
    :param c_config: C Config structure containing configuration parameters.
    :param py_config: Python configuration object from unidecstructure.
    :return: Input structure ready for use in C functions.
    """
    c_input = Input()

    # dataobj.data2 maps to dataMZ and dataInt
    if dataobj.data2 is not None:
        mz = dataobj.data2[:,0]
        intensity = dataobj.data2[:,1]
        c_input.dataMZ = (ctypes.c_float * len(mz))(*mz)
        c_input.dataInt = (ctypes.c_float * len(intensity))(*intensity)
        setattr(c_config, "lengthmz", ctypes.c_int(int(len(mz))))
    else:
        raise ValueError("data2 key not found in dataobj")

    if py_config.masslist is not None:
        # Convert masslist to testmasses
        testmasses = np.array(py_config.masslist, dtype=np.float32)
        c_input.testmasses = (ctypes.c_float * len(testmasses))(*testmasses)
        setattr(c_config, "mfilelen", len(testmasses))
    else:
        # set testmasses to Null if masslist is not provided
        c_input.testmasses = ctypes.POINTER(ctypes.c_float)()
        setattr(c_config, "mfilelen", 0)

    return c_input

"""
struct Decon {
    float *fitdat;
    float *baseline;
    float *noise;
    float *massgrid;
    float *massaxis;
    float *massaxisval;
    float *blur;
    float *newblur;
    float *peakx;
    float *peaky;
    float *dscores;
    int *starttab;
    int *endtab;
    float *mzdist;
    float *rmzdist;
    float error;
    float rsquared;
    int iterations;
    float uniscore;
    float conv;
    float threshold;
    int mlen;
    int plen;
    int scanindex;
};
"""

class Decon(ctypes.Structure):
    _fields_ = [
        ("fitdat", ctypes.POINTER(ctypes.c_float)),
        ("baseline", ctypes.POINTER(ctypes.c_float)),
        ("noise", ctypes.POINTER(ctypes.c_float)),
        ("massgrid", ctypes.POINTER(ctypes.c_float)),
        ("massaxis", ctypes.POINTER(ctypes.c_float)),
        ("massaxisval", ctypes.POINTER(ctypes.c_float)),
        ("blur", ctypes.POINTER(ctypes.c_float)),
        ("newblur", ctypes.POINTER(ctypes.c_float)),
        ("peakx", ctypes.POINTER(ctypes.c_float)),
        ("peaky", ctypes.POINTER(ctypes.c_float)),
        ("dscores", ctypes.POINTER(ctypes.c_float)),
        ("starttab", ctypes.POINTER(ctypes.c_int)),
        ("endtab", ctypes.POINTER(ctypes.c_int)),
        ("mzdist", ctypes.POINTER(ctypes.c_float)),
        ("rmzdist", ctypes.POINTER(ctypes.c_float)),
        ("error", ctypes.c_float),
        ("rsquared", ctypes.c_float),
        ("iterations", ctypes.c_int),
        ("uniscore", ctypes.c_float),
        ("conv", ctypes.c_float),
        ("threshold", ctypes.c_float),
        ("mlen", ctypes.c_int),
        ("plen", ctypes.c_int),
        ("scanindex", ctypes.c_int)
    ]

"""
struct Config {
    char infile[500];
    char outfile[500];
    int numit;
    int numz;
    int endz;
    int startz;
    float zsig;
    float psig;
    float beta;
    float mzsig;
    float msig;
    float molig;
    float massub;
    float masslb;
    int psfun;
    int zpsfun;
    float psmzthresh;
    float mtabsig;
    char mfile[500];
    char manualfile[500];
    int mflag;
    float massbins;
    int limitflag;
    float psthresh;
    int speedyflag;
    int linflag;
    int aggressiveflag;
    float adductmass;
    int rawflag;
    float nativezub;
    float nativezlb;
    int poolflag;
    int manualflag;
    float intthresh;
    float peakshapeinflate;
    int fixedmassaxis;
    int isotopemode;
    int filetype;
    int imflag;
    //IM Parameters
    float dtsig;
    float csig;
    float ccsub;
    float ccslb;
    float ccsbins;
    float temp;
    float press;
    float volt;
    float tcal1;
    float tcal2;
    float tcal3;
    float tcal4;
    int twaveflag;
    float hmass;
    float to;
    float len;
    float edc;
    float nativeccsub;
    float nativeccslb;
    int baselineflag;
    int noiseflag;
    int zout;
    int metamode;
    float minmz;
    float maxmz;
    float mzres;
    float mzbins;
    float bsub;
    float datareduction;
    float peakwin;
    float peakthresh;
    float exwindow;
    int exchoice;
    int exchoicez;
    float exthresh;
    int exnorm;
    int exnormz;
    int peaknorm;
    int orbimode;
    int datanorm;
    //Experimental Parameters
    int filterwidth;
    float zerolog;
    int lengthmz;
    int mfilelen;
    int isolength;
    // DoubleDec Parameters
    int doubledec;
    char kernel[500];
    hid_t file_id;
    char dataset[1024];
    int silent;
    int cdmsflag;
};
"""

class Config(ctypes.Structure):
    _fields_ = [
        ("infile", ctypes.c_char * 500),
        ("outfile", ctypes.c_char * 500),
        ("numit", ctypes.c_int),
        ("numz", ctypes.c_int),
        ("endz", ctypes.c_int),
        ("startz", ctypes.c_int),
        ("zsig", ctypes.c_float),
        ("psig", ctypes.c_float),
        ("beta", ctypes.c_float),
        ("mzsig", ctypes.c_float),
        ("msig", ctypes.c_float),
        ("molig", ctypes.c_float),
        ("massub", ctypes.c_float),
        ("masslb", ctypes.c_float),
        ("psfun", ctypes.c_int),
        ("zpsfun", ctypes.c_int),
        ("psmzthresh", ctypes.c_float),
        ("mtabsig", ctypes.c_float),
        ("mfile", ctypes.c_char * 500),
        ("manualfile", ctypes.c_char * 500),
        ("mflag", ctypes.c_int),
        ("massbins", ctypes.c_float),
        ("limitflag", ctypes.c_int),
        ("psthresh", ctypes.c_float),
        ("speedyflag", ctypes.c_int),
        ("linflag", ctypes.c_int),
        ("aggressiveflag", ctypes.c_int),
        ("adductmass", ctypes.c_float),
        ("rawflag", ctypes.c_int),
        ("nativezub", ctypes.c_float),
        ("nativezlb", ctypes.c_float),
        ("poolflag", ctypes.c_int),
        ("manualflag", ctypes.c_int),
        ("intthresh", ctypes.c_float),
        ("peakshapeinflate", ctypes.c_float),
        ("fixedmassaxis", ctypes.c_int),
        ("isotopemode", ctypes.c_int),
        ("filetype", ctypes.c_int),
        ("imflag", ctypes.c_int),
        # IM Parameters
        ("dtsig", ctypes.c_float),
        ("csig", ctypes.c_float),
        ("ccsub", ctypes.c_float),
        ("ccslb", ctypes.c_float),
        ("ccsbins", ctypes.c_float),
        ("temp", ctypes.c_float),
        ("press", ctypes.c_float),
        ("volt", ctypes.c_float),
        ("tcal1", ctypes.c_float),
        ("tcal2", ctypes.c_float),
        ("tcal3", ctypes.c_float),
        ("tcal4", ctypes.c_float),
        ("twaveflag", ctypes.c_int),
        ("hmass", ctypes.c_float),
        ("to", ctypes.c_float),
        ("len", ctypes.c_float),
        ("edc", ctypes.c_float),
        ("nativeccsub", ctypes.c_float),
        ("nativeccslb", ctypes.c_float),
        ("baselineflag", ctypes.c_int),
        ("noiseflag", ctypes.c_int),
        ("zout", ctypes.c_int),
        ("metamode", ctypes.c_int),
        ("minmz", ctypes.c_float),
        ("maxmz", ctypes.c_float),
        ("mzres", ctypes.c_float),
        ("mzbins", ctypes.c_float),
        ("bsub", ctypes.c_float),
        ("datareduction", ctypes.c_float),
        ("peakwin", ctypes.c_float),
        ("peakthresh", ctypes.c_float),
        ("exwindow", ctypes.c_float),
        ("exchoice", ctypes.c_int),
        ("exchoicez", ctypes.c_int),
        ("exthresh", ctypes.c_float),
        ("exnorm", ctypes.c_int),
        ("exnormz", ctypes.c_int),
        ("peaknorm", ctypes.c_int),
        ("orbimode", ctypes.c_int),
        ("datanorm", ctypes.c_int),
        # Experimental Parameters
        ("filterwidth", ctypes.c_int),
        ("zerolog", ctypes.c_float),
        ("lengthmz", ctypes.c_int),
        ("mfilelen", ctypes.c_int),
        ("isolength", ctypes.c_int),
        # DoubleDec Parameters
        ("doubledec", ctypes.c_int),
        ("kernel", ctypes.c_char * 500),
        ("file_id", ctypes.c_int),  # hid_t, adjust if needed
        ("dataset", ctypes.c_char * 1024),
        ("silent", ctypes.c_int),
        ("cdmsflag", ctypes.c_int),
    ]


unideclib = ctypes.CDLL(default_dll_path)
unideclib.run_unidec_core.argtypes = [Config, Input, ctypes.POINTER(Decon), ctypes.c_int, ctypes.c_int]
unideclib.run_unidec_core.restype = ctypes.c_int

unideclib.SetupZtab.argtypes = [Config, ctypes.POINTER(Input)]
unideclib.SetupZtab.restype = ctypes.c_void_p

unideclib.PostImport.argtypes = [ctypes.POINTER(Config)]
unideclib.PostImport.restype = ctypes.c_void_p

unideclib.SetDefaultConfig.argtypes = [ctypes.POINTER(Config)]
unideclib.SetDefaultConfig.restype = ctypes.c_void_p

def py_to_c_config(py_config):
    """
    Convert a Python config dictionary to a C Config structure.
    :param py_config: Dictionary containing configuration parameters.
    :return: Config structure ready for use in C functions.
    """
    c_config = Config()
    unideclib.SetDefaultConfig(ctypes.byref(c_config))
    for field in c_config._fields_:
        name, ctype = field
        if name in py_config:
            value = py_config[name]
            value, type = value
            if type == str:
                value = str(value).encode('utf-8')
            elif type == float:
                value = ctypes.c_float(float(value))
            elif type == int:
                value = ctypes.c_int(int(value))
            try:
                setattr(c_config, name, value)
            except TypeError as e:
                print(f"Error setting field {name} with value {value}. Expected type {ctype}.")
                raise e

    return c_config


def c_to_py_decon(c_decon, eng):
    """
    Convert a C Decon structure to a Python dictionary.
    :param c_decon: C Decon structure containing deconvolution results.
    :param eng: UniDec engine instance for context.
    :return: Dictionary containing deconvolution results.
    """
    mlen = c_decon.mlen
    mzlen = len(eng.data.data2[:, 0]) if eng.data.data2 is not None else 0
    zlen = eng.config.numz
    eng.data.fitdat = np.ctypeslib.as_array(c_decon.fitdat, shape=(mzlen,))
    massx = np.ctypeslib.as_array(c_decon.massaxis, shape=(mlen,))
    massy = np.ctypeslib.as_array(c_decon.massaxisval, shape=(mlen,))
    eng.data.massdat = np.column_stack((massx, massy))

    eng.data.massgrid = np.ctypeslib.as_array(c_decon.massgrid, shape=(mlen * zlen,))

    eng.data.mzgrid = np.ctypeslib.as_array(c_decon.blur, shape=(mzlen * zlen,))


def run_unidec_core(eng):
    starttime = time.perf_counter()
    c_config = py_to_c_config(eng.config)
    unideclib.PostImport(ctypes.byref(c_config))
    c_input = py_to_c_input(eng.data, c_config, eng.config)
    c_decon = Decon()
    # Setup Ztab
    unideclib.SetupZtab(c_config, ctypes.byref(c_input))
    # Run the core UniDec function
    result = unideclib.run_unidec_core(c_config, c_input, ctypes.byref(c_decon), 0, 0)
    if result != 0:
        raise RuntimeError(f"UniDec core function failed with error code {result}")
    c_to_py_decon(c_decon, eng)
    print(f"UniDec core processing time: {time.perf_counter() - starttime:.2f} seconds")


if __name__ == "__main__":
    eng = UniDec()
    testfile = os.path.join(os.path.join(os.path.dirname(default_dll_path), "Example Data"), "ADH.txt")
    eng.open_file(testfile)

    run_unidec_core(eng)

    # c_config = py_to_c_config(eng.config)
    # c_input = py_to_c_input(eng.data)
    # wrapper.close()
