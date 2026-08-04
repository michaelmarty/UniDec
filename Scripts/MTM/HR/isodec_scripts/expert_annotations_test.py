import torch
from rdkit.sping.colors import transparent

from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection, compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pyperclip
from unidec.IsoDec.match import *
from Scripts.MTM.proteomicstools import *
from unidec.UniDecImporter.Importer import ImporterFactory

torch.cuda.is_available()

# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\PXD018043\\aTpm_ECD_20170616_F1-8_042917_ECD_1364.4mz_3width_0.7V_20ms_1.0s_1000scans_000001.mzML"
os.chdir(os.path.dirname(file2))

mzrange = [500, 600]
mzrange = None


eng = IsoDecEngine(phaseres=8, verbose=False, use_wrapper=False)

eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

reader = ImporterFactory.create_importer(file2)
# spectrum = reader.grab_centroid_data(1)
spectrum = reader.grab_scan_data(0)
# spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

    print("Spectrum Length:", len(spectrum))

expert_peaks = read_msalign_to_matchedcollection(
    "Z:\\Group Share\\JGP\\PXD018043\\expert_annotations_2020_04_29\\expert_annotations_2020_04_29\\aTpm\\ECD\\Expert Annotations\\aTpm_ECD_Ground_Truth.msalign",
    spectrum)
if mzrange is not None:
    expert_peaks = expert_peaks.filter(mzrange[0], mzrange[1], type="mz")

expert_peaks.merge_missed_monoisotopics()

timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True, window=80)

pyperclip.copy(eng.pks.copy_to_string())

shared, u1, u2 = compare_matchedcollections(eng.pks, expert_peaks, ppmtol=50)

print("Time:", time.time() - timestart)

print("N Peaks:", len(eng.pks.peaks))
print("Shared Peaks: " + str(len(shared)), "Unique Peaks IsoDec: " + str(len(u1)),
      "Unique Peaks Manual: " + str(len(u2)))

try:
    percent_matched = len(shared) / len(expert_peaks.peaks) * 100
    print("Percent of annotated peaks matched: " + str(percent_matched))

    percent_clusters_matched = len(shared) / len(eng.pks.peaks) * 100
    print("Percent of accepted clusters with a match: " + str(percent_clusters_matched))
except:
    pass
# set axis fonts to 12
tickfont = 12
labfont = 12
plot_pks(u1, forcecolor="r", show=False, plotmass=False, nocentroids=True, tickfont=tickfont, labfont=labfont)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True, tickfont=tickfont, labfont=labfont)
plot_pks(shared, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=True,
         tickfont=tickfont, labfont=labfont)

# outdir = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures"
# os.chdir(outdir)
# plt.savefig("expert_annotations.png", dpi=300, transparent=True)

plt.show()
