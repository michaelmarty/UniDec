import ctypes
import os
import numpy as np
import time
import platform
from unidec.tools import start_at_iso

current_path = os.path.dirname(os.path.realpath(__file__))

if platform.system() == "Windows":
    dllname = "isogenmass.dll"
    fftdllname = "isofft.dll"
elif platform.system() == "Linux":
    dllname = "isogenmass.so"
    fftdllname = "isofft.so"
else:
    dllname = "isogenmass.dylib"
    fftdllname = "isofft.dylib"



dllpath = start_at_iso(dllname, guess=current_path)
if not dllpath:
    print("DLL not found anywhere:", dllname)
# else:
#     print("DLL Found at:", dllpath)

isodist = ctypes.c_float * 64

isogen_c_lib = ctypes.CDLL(dllpath)
isogen_c_lib.isogenmass.argtypes = [
            ctypes.c_float,
            ctypes.POINTER(ctypes.c_float)
        ]

isogenmass = isogen_c_lib.isogenmass

isogen_c_lib.isomike.argtypes = [
            ctypes.c_float,
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_int,
        ]

isomike = isogen_c_lib.isomike

dllpath2 = start_at_iso(fftdllname, guess=current_path)
# dllpath2 = "C:\\Python\\UniDec3\\Scripts\\MTM\\IsoFFT\\isofft.dll"
if not os.path.isfile(dllpath2):
    raise Exception("DLL not found")
isofft_c_lib = ctypes.CDLL(dllpath2)

isofft = isofft_c_lib.isodist_from_averagine_mass
isofft.argtypes = [
    ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.c_int,
]


def gen_isofft(mass, isolen=64):
    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    # Call the C function
    isofft(
        ctypes.c_float(mass),
        isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        ctypes.c_int(isolen),
        ctypes.c_int(0),
    )
    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist

def gen_isodist(mass, isolen=64):
    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    # Call the C function
    isogenmass(
        ctypes.c_float(mass),
        isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    )
    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist

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
    def gen_isodist(self, mass):
        return gen_isodist(mass, self.isolen)

    def gen_isomike(self, mass):
        return gen_isomike(mass, self.isolen)

    def gen_isofft(self, mass):
        return gen_isofft(mass, self.isolen)

if __name__ == "__main__":

    m =20000
    d1 = gen_isodist(m)
    d2 = gen_isomike(m)
    d3 = gen_isofft(m)

    import matplotlib.pyplot as plt
    plt.plot(d1/np.amax(d1), label="Isogen")
    plt.plot(d2/np.amax(d2), label="Isomike")
    plt.plot(d3/np.amax(d3), label="FFT")
    plt.legend()
    plt.show()


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
