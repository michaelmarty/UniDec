import os
from unidec.IsoDec.engine import IsoDecEngine
import unidec.tools as ud
import time

tfd = 21.55

tstart = time.perf_counter()

topdir = "Z:\\Group Share\\JGP\\MSV000084001_FlashDeconv\\raw"

os.chdir(topdir)

rawfiles = ud.match_files(topdir, "*.raw")

eng = IsoDecEngine(phaseres=4, use_wrapper=False)

for file in rawfiles:
    print(file)
    eng.process_file(file, scans=None)
    eng.export_peaks(type="msalign", filename=file)

print("Time:", time.perf_counter() - tstart)


