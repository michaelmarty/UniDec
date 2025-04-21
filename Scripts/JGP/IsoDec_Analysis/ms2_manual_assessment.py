import os
import matplotlib as mpl

from Scripts.MTM.proteomicstools import searchTD
from unidec.IsoDec.engine import IsoDecEngine
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.IsoDec.match import *
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *

valid_count = 0
invalid_count = 0

valid_region = (0,10000,0,1000000)
invalid_region = (0,10000,-1000000,0)

def on_click(event):
    global valid_count, invalid_count
    x = event.xdata
    y = event.ydata

    if x is None or y is None:
        return

    if valid_region[0] <= x <= valid_region[1] and valid_region[2] <= y <= valid_region[3]:
        valid_count += 1
        print("Valid peak, count:", valid_count)
    elif invalid_region[0] <= x <= invalid_region[1] and invalid_region[2] <= y <= invalid_region[3]:
        invalid_count += 1
        print("Invalid peak, count:", invalid_count)

def perform_TD_assessment(id8_mc, fd_mc, tfd_mc, pfm_sequence, spectrum):
    print("Performing TD Assessment")
    searchTD(id8_mc.peaks, pfm_sequence, ["b", "y", "c","z"], "ppm", 20, verbose=True)
    searchTD(fd_mc.peaks, pfm_sequence, ["b", "y", "c","z"], "ppm", 20, verbose=True)
    searchTD(tfd_mc.peaks, pfm_sequence, ["b", "y", "c","z"], "ppm", 20, verbose=True)

    #Remove all peaks with no matchedions
    id8_mc.peaks = [p for p in id8_mc.peaks if len(p.matchedions) > 0]
    fd_mc.peaks = [p for p in fd_mc.peaks if len(p.matchedions) >0]
    tfd_mc.peaks = [p for  p  in tfd_mc.peaks if len(p.matchedions) > 0]

    shared_id8_FD, uid8_fd, ufd_id8 = compare_matchedcollections(id8_mc, fd_mc, ppmtol=50, maxshift=3)
    shared_id8_tfd, uid8_tfd, utfd_id8 = compare_matchedcollections(id8_mc, tfd_mc, ppmtol=50, maxshift=3)

    uid8, u1, u2 = compare_matchedcollections(uid8_fd, uid8_tfd, ppmtol=50, maxshift=3)

    plot_pks(uid8, forcecolor="r", centroids=spectrum, ccolor="k", show=False, plotmass=False, nocentroids=True)
    plt.show()



#Need to load the .raw file
#Need to load the .msalign file for each algorithm
#Need to plot the peaks for the selected scan
#Need a way to manually validate or invalidate the peaks

# Set backend to Agg
mpl.use('WxAgg')

rawfile = r"Z:\Group Share\JGP\PXD019247\20170105_L_MaD_ColC4_Ecoli20161108-11_01.raw"

id8_msalign = r"Z:\Group Share\JGP\PXD019247\TD_Search\msalign_for_toppic\8bin\20170105_L_MaD_ColC4_Ecoli20161108-11_01.msalign_ms2.msalign"
fd_msalign = r"Z:\Group Share\JGP\PXD019247\TD_Search\FLASHDeconv_msalign\20170105_L_MaD_ColC4_Ecoli20161108-11_01_ms2.msalign"
tfd_msalign = r"Z:\Group Share\JGP\PXD019247\TD_Search\topFD_msalign\20170105_L_MaD_ColC4_Ecoli20161108-11_01_ms2.msalign"

current_scan = 5113

importer = ImporterFactory.create_importer(rawfile)
spectrum = importer.grab_centroid_data(current_scan)

for i in range(len(spectrum)):
    plt.plot([spectrum[i][0], spectrum[i][0]], [0, spectrum[i][1]], color="k")
plt.axhline(0, color="k")
plt.show()

print("Reading msAlign files...")
id8_mc = read_msalign_to_matchedcollection_fast(id8_msalign, data=spectrum)
fd_mc = read_msalign_to_matchedcollection_fast(fd_msalign, data=spectrum)
tfd_mc = read_msalign_to_matchedcollection_fast(tfd_msalign, data=spectrum)

#Remove all peaks in the matched collection that are not in the current scan
id8_mc.peaks = [p for p in id8_mc.peaks if p.scan == current_scan]
fd_mc.peaks = [p for p in fd_mc.peaks if p.scan == current_scan]
tfd_mc.peaks = [p for p in tfd_mc.peaks if p.scan == current_scan]

print("ID8 Peaks:", len(id8_mc.peaks))
print("FD Peaks:", len(fd_mc.peaks))
print("TFD Peaks:", len(tfd_mc.peaks))

pfm_sequence = 'MKAKELREKSVEELNTELLNLLREQFNLRMQAASGQLQQSHLLKQVRRDVARVKTLLNEKAGA'

perform_TD_assessment(id8_mc, fd_mc, tfd_mc, pfm_sequence, spectrum)
exit()

for p in id8_mc.peaks:
    min_mz = p.mz - 1.05
    max_mz = p.mz + 2.05

    if min_mz > np.min(p.isodist[:,0] - 1.5):
        min_mz = np.min(p.isodist[:,0] - 1.5)
    if max_mz < np.max(p.isodist[:,0]) + 1.5:
        max_mz = np.max(p.isodist[:,0] + 1.5)

    #Get the spectrum in the range of the peak
    spec = ud.datachop(spectrum, min_mz, max_mz)
    for i in range(len(spec)):
        plt.plot([spec[i][0], spec[i][0]], [0, spec[i][1]], color="k")
    for i in range(len(p.isodist)):
        plt.plot([p.isodist[i][0], p.isodist[i][0]], [0, -p.isodist[i][1]], color="r")

    #Plot a horizontal line at 0
    plt.axhline(0, color="k")
    cid = plt.gcf().canvas.mpl_connect('button_press_event', on_click)
    plt.show()


print("Valid peaks:", valid_count)
print("Invalid peaks:", invalid_count)