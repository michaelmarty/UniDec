import os
import fnmatch
import numpy as np
from unidec.IsoDec.runtime import IsoDecRuntime

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

filepath = r"C:\Users\Admin\Desktop\martylab\IsoDec\PXD019247\PXD019247\20170105_L_MaD_ColC4_Ecoli20161108-10_01.raw"
os.chdir(os.path.dirname(filepath))

eng = IsoDecRuntime(phaseres=8)

topdir = r"C:\Users\Admin\Desktop\martylab\IsoDec\PXD019247\PXD019247"
rawfiles = match_files(topdir, "*.raw")
print("N Raw Files:", len(rawfiles))
for file in rawfiles:
    print(file)

for file in rawfiles:
    filename = os.path.splitext(file)[0]
    print("Processing:", filename)
    eng.process_file(file)

    eng.export_peaks(type="msalign", filename=filename, reader=eng.reader, act_type="EThcD", max_precursors=1)