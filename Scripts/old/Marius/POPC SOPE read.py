import numpy
import numpy as np

import unidec.tools as ud
from unidec.modules import peakstructure


#set the average lipid mass and sub
def sub_div(eng, sub=44088, div=756.5):
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

actualdata = numpy.loadtxt(r'Z:\Group Share\Marius Kostelic\EVEN More mixed lipids\figures\figure 2\mammalian\POPC SOPE.csv', delimiter = ",")
print(actualdata[1])
mylist=numpy.mean[actualdata[:3]]
print(mylist)

