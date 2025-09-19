import os
import numpy
import numpy as np

import unidec.tools as ud
from unidec.metaunidec.mudeng import MetaUniDec
from unidec.modules import peakstructure


#set the average lipid mass and sub
def sub_div(eng, sub=44088, div=726.0):
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
    return sd_results



#set the path or directory
path = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\Bacterial\\DPPG SOPE POCL"
os.chdir(path)
#read HDF5 files in folder
files = ud.match_files(path, "*.hdf5")
print(files)

lipids = ["DPPG_SOPE_POCL"]
modes = ["pos", "neg", "imid"]#"methyl"]
replicates = ["1", "2", "3"]
outdata = list()
eng = MetaUniDec()
#Runs through each file and prints average lipid masses, lipid name, and replicate number
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
                sd_results = sub_div(eng)
                outdata.append(sd_results)
                #outdata.append(l)
                #outdata.append(m)
                #outdata.append(r)
                print(l,m,r)
#for outstring in replicates

#a = numpy.asarray([outdata])
#saves data as a csv file
numpy.savetxt(r'Z:\Group Share\Marius Kostelic\EVEN More mixed lipids\Bacterial\DPPG SOPE POCL\DPPG SOPE POCL all.csv', outdata, delimiter = ",")
#actualdata = numpy.loadtxt(r'Z:\Group Share\Marius Kostelic\EVEN More mixed lipids\figures\figure 2\mammalian\POPC SOPE.csv', delimiter = ",")
#print(actualdata[1])

#files = ud.match_files(path, "*_.hdf5")
#filename = "20191030_MK_POPC_SOPE_2methylImid_Ramp_0_200_20.hdf5"
#df = h5py.File(filename, "r+")
#print(df)



