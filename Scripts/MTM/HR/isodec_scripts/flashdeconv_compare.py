from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection, compare_matchedcollections
import time
from unidec.IsoDec.match import *

timing = {"FlashDeconv": 1295.0, "TopFD": 637, "IsoDec Python": 298, "IsoDec C": 118.0}

tstart = time.perf_counter()

topdir = "Z:\\Group Share\\JGP\\MSV000084001_FlashDeconv\\raw"

fddir = "Z:\\Group Share\\JGP\\MSV000084001_FlashDeconv\\output"

os.chdir(topdir)

rawfiles = ud.match_files(topdir, "*.raw")


for file in rawfiles:
    print(file)
    isodecfile = file.replace(".raw", "_ms1.msalign")
    print(isodecfile)
    idpeaks = read_msalign_to_matchedcollection(isodecfile)
    print(len(idpeaks))


    fdfile = os.path.join(fddir, os.path.basename(file).replace(".raw", "_ms1.msalign"))
    print(fdfile)
    fdpeaks = read_msalign_to_matchedcollection(fdfile)
    print(len(fdpeaks))

    shared, u1, u2 = compare_matchedcollections(idpeaks, fdpeaks, ppmtol=10)

    print("Shared Peaks: " + str(len(shared)), "Unique Peaks IsoDec: " + str(len(u1)),
          "Unique Peaks FlashDeconv: " + str(len(u2)))

    try:
        percent_matched = (len(fdpeaks) - len(u2)) / len(fdpeaks) * 100
        print("Percent of fd peaks matched: " + str(percent_matched))

        percent_fd_unique = len(u2) / len(fdpeaks) * 100
        print("Percent of fd peaks unique: " + str(percent_fd_unique))

        percent_id_matched = (len(idpeaks) - len(u1)) / len(idpeaks) * 100
        print("Percent of id peaks matched: " + str(percent_id_matched))

        percent_id_unique = len(u1) / len(idpeaks) * 100
        print("Percent of id peaks unique: " + str(percent_id_unique))
    except:
        pass

print("Time:", time.perf_counter() - tstart)


