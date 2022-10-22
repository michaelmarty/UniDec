from pyopenms import *
import os
import unidec_modules.unidectools as ud
from copy import deepcopy


def fix_parent_mz(exp, i, p, window=2, purity_thresh=0.5):
    j = deepcopy(i)
    isMS1 = False
    while not isMS1:
        if exp[j].getMSLevel() == 1:
            isMS1 = True
        else:
            j -= 1
    ms1 = exp[j]

    pmass = p.getMZ()
    mz, intensity = ms1.get_peaks()
    data = np.transpose([mz, intensity])
    chopdat = ud.datachop(data, pmass - window, pmass + window)
    if len(chopdat) > 0:
        index = np.argmax(chopdat[:, 1])
        peakmz = chopdat[index, 0]
        peakint = chopdat[index, 1]
        print(pmass, peakmz, peakint)

        p.setMZ(peakmz)
        p.setIntensity(peakint)
        p.setCharge(1)

        p.setIsolationWindowLowerOffset(window)
        p.setIsolationWindowUpperOffset(window)

        if purity_thresh is not None:
            purity_score = PrecursorPurity().computePrecursorPurity(ms1, p, 10, True)
            prop = purity_score.signal_proportion
            print("Proportion:", prop)
            if prop < purity_thresh:
                return None

        polarity = ms1.getInstrumentSettings().getPolarity()
        if polarity == IonSource.Polarity.POSITIVE:
            p.setCharge(1)
        elif polarity == IonSource.Polarity.NEGATIVE:
            p.setCharge(-1)

        return p
    else:
        return None


toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

mzML_files = ud.match_files(toppath, "*.mzml", exclude="fixed")

for file in mzML_files:
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    filtered = MSExperiment()
    for i, s in enumerate(exp):
        if s.getMSLevel() == 2:
            p = s.getPrecursors()[0]
            p = fix_parent_mz(exp, i, p)
            if p is not None:
                s.setPrecursors([p])
                filtered.addSpectrum(s)
        else:
            filtered.addSpectrum(s)

    MzMLFile().store(file[:-5] + "_fixed.mzML", filtered)
exit()
