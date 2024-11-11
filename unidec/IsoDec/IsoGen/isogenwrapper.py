import ctypes
import os
import numpy as np
import time

from unidec.tools import start_at_iso

dllpath = start_at_iso("isogenmass.dll")
if not dllpath:
    print("DLL not found anywhere")


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
        # Create empty array
        isodist = np.zeros(self.isolen).astype(np.float32)
        # Call the C function
        isogenmass(
            ctypes.c_float(mass),
            isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        )
        # Convert isodist to numpy
        isodist = np.ctypeslib.as_array(isodist)
        return isodist

    def gen_isomike(self, mass):
        # Create empty array
        isodist = np.zeros(self.isolen).astype(np.float32)
        # Call the C function
        isomike(
            ctypes.c_float(mass),
            isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(self.isolen),
            ctypes.c_int(0),
        )
        # Convert isodist to numpy
        isodist = np.ctypeslib.as_array(isodist)
        return isodist

# if __name__ == "__main__":
#
#     m =1000
#     gen_isodist(m)
#     gen_isomike(m)
#     exit()
#
#     eng = IsoGenWrapper(dllpath=dllpath)
#     n = 1000
#     random_masses = np.random.uniform(1000, 10000, n)
#     starttime = time.perf_counter()
#     for mass in random_masses:
#         dist = eng.gen_isomike(mass)
#
#     print("Time:", time.perf_counter()-starttime)
#     print("Time Per:", (time.perf_counter()-starttime)/n)
