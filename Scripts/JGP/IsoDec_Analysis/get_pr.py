from unidec.IsoDec.compare_tools import read_manual_annotations, compare_annotated
from unidec.IsoDec.engine import *

# Set backend to Agg
mpl.use('WxAgg')

etd_spectrum = np.loadtxt(r"Z:\Group Share\JGP\js8b05641_si_001\etd_spectrum.txt")
uvpd_spectrum = np.loadtxt(r"Z:\Group Share\JGP\JPST001885\UVPD_Test\uvpd_spectrum.txt")
tof_spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt")

etd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt", delimiter=' ')
uvpd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt", delimiter='\t')
tof_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt", delimiter='\t')

total_annotated = len(etd_annotated) + len(uvpd_annotated) + len(tof_annotated)
print("Total Annotations:", total_annotated)

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
# eng.config.css_thresh = 0.75
# eng.config.minareacovered = 0.2
# eng.config.knockdown_rounds = 4

truepos = 0
falsepos = 0
falsneg = 0

eng.batch_process_spectrum(etd_spectrum, centroided=True)
shared, u1, u2 = compare_annotated(etd_annotated, eng.pks, ppmtol=50, maxshift=3)
truepos += len(shared)
falsepos += len(u2)
falsneg += len(u1)

eng.pks = MatchedCollection()

eng.batch_process_spectrum(uvpd_spectrum, centroided=True)
shared, u1, u2 = compare_annotated(uvpd_annotated, eng.pks, ppmtol=50, maxshift=3)
truepos += len(shared)
falsepos += len(u2)
falsneg += len(u1)

# eng.pks = MatchedCollection()
#
# eng.batch_process_spectrum(tof_spectrum, centroided=True)
# shared, u1, u2 = compare_annotated(tof_annotated, eng.pks, ppmtol=50, maxshift=3)
# truepos += len(shared)
# falsepos += len(u2)
# falsneg += len(u1)

print("True Positives:", truepos)
print("False Positives:", falsepos)
print("False Negatives:", falsneg)
print("Precision:", truepos/(truepos+falsepos))
print("Recall:", truepos/(truepos+falsneg))
