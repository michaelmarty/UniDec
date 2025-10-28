import numpy as np
import os
import unidec.tools as ud
from unidec.metaunidec.mudeng import MetaUniDec
from unidec.modules import peakstructure


#set the average lipid mass and sub
def sub_div(eng, sub=44088, div=760):
    sd_results = []
    message = "Average Masses:\n\n"
    outstring = ""
    for s in eng.data.spectra:
        pks = peakstructure.Peaks()
        peaks = ud.peakdetect(s.massdat, eng.config)
        pks.add_peaks(peaks, massbins=eng.config.massbins)
        sd_result = ud.subtract_and_divide(pks, sub, div)
        sd_results.append(sd_result)

        outstring += str(sd_result) + "\n"

    avg = np.mean(sd_results)
    message += outstring
    message += "\nOverall Average: " + str(avg)
    outstring += "\n" + str(avg)
    print(message, outstring)

#set the path or directory
path = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\MAMMALIAN\\POPC"
os.chdir(path)
#read HDF5 files in folder
files = ud.match_files(path, "*.hdf5")
print(files)

lipids=["POPC"]
modes=["pos", "neg", "imid", "methyl"]
replicates = ["1", "2", "3"]

eng = MetaUniDec()

for l in lipids:
    for m in modes[:3]:
        for r in replicates[:3]:
            filename=None
            for f in files:
                if l in f:
                    if m in f:
                        if str(r) in f:
                            filename=f
                            print(f)
            if filename is not None:
                eng.open(filename)
                sub_div(eng)
                print(l,m,r)
#for outstring in replicates


#files = ud.match_files(path, "*_.hdf5")
#filename = "20191030_MK_POPC_SOPE_2methylImid_Ramp_0_200_20.hdf5"
#df = h5py.File(filename, "r+")
#print(df)

