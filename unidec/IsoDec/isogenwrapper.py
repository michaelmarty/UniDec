import ctypes
import os
import numpy as np
import time
import platform
from unidec.tools import start_at_iso

current_path = os.path.dirname(os.path.realpath(__file__))


#Replace all of this with the isogen.dll for both isogenmass AND isogen
if platform.system() == "Windows":
    dllname = "isogen.dll"
elif platform.system() == "Linux":
    dllname = "libisogen.so"

else:
    print("Not yet implemented for MacOS")
    dllname = "isogen.dylib"



#need to commmit isogendll to srccmake
dllpath = start_at_iso(dllname, guess=current_path)

if not dllpath:
    print("DLL not found anywhere:", dllname)
# else:
#     print("DLL Found at:", dllpath)

isodist = ctypes.c_float * 64


isogen_c_lib = ctypes.CDLL(dllpath)
# isogen_c_lib.isogenmass.argtypes = [
#             ctypes.c_float,
#             ctypes.POINTER(ctypes.c_float)
#         ]

isogen_c_lib.fft_rna_mass_to_dist.argtypes = [ctypes.c_float]
isogen_c_lib.fft_rna_mass_to_dist.restype = ctypes.POINTER(ctypes.c_float)
isogen_c_lib.fft_pep_mass_to_dist.restype = ctypes.c_float
isogen_c_lib.fft_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]


isogen_c_lib.isogen_nn.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_char_p]
isogen_c_lib.isogen_nn.restype = None

isogen_c_lib.isomike.argtypes = [
            ctypes.c_float,
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_int,
        ]

isomike = isogen_c_lib.isomike


def gen_isodist_nn(mass, type, isolen=64):
    res = b""
    if type is None:
        return
    if type == "RNA":
        res = b"Rna"
    elif type == "PEPTIDE":
        res = b"Peptide"
        # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    # Call the C function
    isogen_c_lib.isogen_nn(
        ctypes.c_float(mass),
        ptr,
        ctypes.c_int(isolen),
        res,
    )
    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist

def fft_gen_isodist(mass, type, isolen = 128):
    if type == "RNA":
        mass_c = ctypes.c_float(mass)
        isolist_ptr = isogen_c_lib.fft_rna_mass_to_dist(mass_c)
        isolist = np.ctypeslib.as_array(isolist_ptr, shape=(128,))
        return isolist
    elif type == "PEPTIDE":
        fmass = ctypes.c_float(mass)
        isodist = (ctypes.c_float * isolen)()
        offset = ctypes.c_int(0)
        isolen = ctypes.c_int(isolen)
        isogen_c_lib.fft_pep_mass_to_dist(fmass, isodist, isolen, offset)
        return np.array(isodist).copy()



def gen_isomike(mass, isolen=64):
    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    # Call the C function
    isomike(
        ctypes.c_float(mass),
        isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        ctypes.c_int(isolen),
        ctypes.c_int(0),
    )
    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist


class IsoGenWrapper:
    def __init__(self, dllpath=None):
        self.isolen = 64
    def gen_isodist(self, mass, type="PEPTIDE"):
        return gen_isodist_nn(mass, type, self.isolen)

    def gen_isomike(self, mass):
        return gen_isomike(mass, self.isolen)

    def gen_isodist_fft(self, mass):
        return fft_gen_isodist(mass, type, self.isolen)

if __name__ == "__main__":

    m =20000
    d1 = gen_isodist(m)
    d2 = gen_isomike(m)

    import matplotlib.pyplot as plt
    plt.plot(d1/np.amax(d1), label="Isogen")
    plt.plot(d2/np.amax(d2), label="Isomike")
    plt.legend()
    #plt.show()


    eng = IsoGenWrapper(dllpath=dllpath)
    n = 10000
    random_masses = np.random.uniform(1000, 60000, n)
    starttime = time.perf_counter()
    for mass in random_masses:
        dist = eng.gen_isodist(mass)

    print("Isogen Time:", time.perf_counter()-starttime)
    print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)

    starttime = time.perf_counter()
    for mass in random_masses:
        dist = eng.gen_isomike(mass)

    print("Isomike Time:", time.perf_counter()-starttime)
    print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)

    starttime = time.perf_counter()
    for mass in random_masses:
        dist = eng.gen_isofft(mass)

    print("FFT Time:", time.perf_counter()-starttime)
    print("Microseconds Per:", (time.perf_counter()-starttime)/n * 1e6)
