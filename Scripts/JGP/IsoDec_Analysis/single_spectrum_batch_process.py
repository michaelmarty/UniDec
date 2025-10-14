import matplotlib.pyplot as plt

from Scripts.MTM.proteomicstools import searchTD
from unidec.IsoDec.compare_tools import read_manual_annotations, compare_annotated, compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.match import *
import bisect
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

class PeakGroup:
    def __init__(self, annotated, peak, source):
        self.annotated = annotated
        self.mz = peak.mz
        self.peaks = np.array([peak])
        self.sources = np.array([source])

def plot_peak_group(spectrum, group):
    mz_range = [group.annotated.mz - 1.5, group.annotated.mz + 1.5]
    chopped = ud.datachop(spectrum, mz_range[0], mz_range[1])
    for i in range(len(chopped)):
        plt.vlines(chopped[i,0],0, chopped[i,1], color='k' )
    for i in range(len(group.annotated.isodist)):
        plt.vlines(group.annotated.isodist[i,0],0, group.annotated.isodist[i,1], color='g')

    for i in range(len(group.peaks[0].isodist)):
        plt.vlines(group.peaks[0].isodist[i,0], -group.peaks[0].isodist[i,1], 0, color='r')


    for i in range(len(group.peaks[1].isodist)):
        plt.vlines(group.peaks[1].isodist[i,0],-group.peaks[1].isodist[i,1], 0, color='b')

    print("Annotated mz: ", group.annotated.mz, "Annotated z: ", group.annotated.z)
    print("Peak 1 mz: ", group.peaks[0].mz, "Peak 1 z: ", group.peaks[0].z, "Peak 1 source: ", group.sources[0])
    print("Peak 2 mz: ", group.peaks[1].mz, "Peak 2 z: ", group.peaks[1].z, "Peak 2 source: ", group.sources[1])
    plt.hlines(0, mz_range[0], mz_range[1], color='k')
    plt.xlim(mz_range[0], mz_range[1])
    plt.show()

def print_precision_recall(shared, u_annotated, u_experimental, existing_truepos, existing_falsepos, existing_falseneg):
    existing_truepos += shared
    existing_falsepos += u_experimental
    existing_falseneg += u_annotated

    precision = shared / (shared + u_experimental)
    recall = shared / (shared + u_annotated)
    print("Precision: ", precision)
    print("Recall: ", recall)
    f1 = (2*precision*recall)/(precision+recall)
    print("F1: ", f1)
    return existing_truepos, existing_falsepos, existing_falseneg

def get_inconsistincies(spectrum, annotated, peaks_nn, peaks_thrash, maxshift=3, ppmtol=20):
    annotated.peaks = sorted(annotated.peaks, key=lambda p: p.mz)
    peaks_nn.peaks = sorted(peaks_nn.peaks, key=lambda p: p.mz)
    peaks_thrash.peaks = sorted(peaks_thrash.peaks, key=lambda p: p.mz)

    annotated_mzs = np.array([p.mz for p in annotated.peaks])
    nn_mzs = np.array([p.mz for p in peaks_nn.peaks])
    thrash_mzs = np.array([p.mz for p in peaks_thrash.peaks])

    groups = np.array([])
    group_mzs = np.array([])

    for p_nn in peaks_nn.peaks:
        #First determine if there's an annotated cluster at the same m/z
        nearest, ann_within_tol = fastwithin_abstol_withnearest(annotated_mzs, p_nn.mz, 0.2)
        if len(ann_within_tol) > 0:
            group = PeakGroup(annotated.peaks[nearest], p_nn, "NN")
            if len(groups) == 0:
                groups = np.append(groups, group)
                group_mzs = np.append(group_mzs, nearest)
            else:
                group_idx = bisect.bisect_left(group_mzs, p_nn.mz)
                groups = np.insert(groups, group_idx, group)
                group_mzs = np.insert(group_mzs, group_idx, p_nn.mz)

    for p_thrash in peaks_thrash.peaks:
        #First check if there's an existing group at this m/z
        nearest, groups_within_tol = fastwithin_abstol_withnearest(group_mzs, p_thrash.mz, 0.2)
        if len(groups_within_tol) > 0:
            groups[nearest].sources = np.append(groups[nearest].sources, "THRASH")
            groups[nearest].peaks = np.append(groups[nearest].peaks, p_thrash)
        else:
            #Check if there's an annotated peak within the tolerance
            nearest, ann_within_tol = fastwithin_abstol_withnearest(annotated_mzs, p_thrash.mz, 0.2)
            if len(ann_within_tol) > 0:
                group = PeakGroup(annotated.peaks[nearest], p_thrash, "THRASH")
                if len(groups) == 0:
                    groups = np.append(groups, group)
                    group_mzs = np.append(group_mzs, nearest)
                else:
                    group_idx = bisect.bisect_left(group_mzs, p_thrash.mz)
                    groups = np.insert(groups, group_idx, group)
                    group_mzs = np.insert(group_mzs, group_idx, p_thrash.mz)

    both_matched = [group for group in groups if len(group.sources) == 2]
    nn_only = [group for group in groups if len(group.sources) == 1 and "NN" in group.sources]
    thrash_only = [group for group in groups if len(group.sources) == 1 and "THRASH" in group.sources]

    for group in both_matched:
        if group.peaks[0].z == group.peaks[1].z:
            continue
        else:
            plot_peak_group(spectrum, group)

eng_nn = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
eng_nn.config.activescan = 1
eng_nn.config.activescanrt = -1
eng_nn.config.activescanorder = 2
eng_nn.config.knockdown_rounds = 4
eng_nn.config.css_thresh = 0.7451
eng_nn.config.minareacovered = 0.15
eng_nn.predmode = 0

nn_true_pos = 0
nn_false_pos = 0
nn_false_neg = 0

eng_thrash = IsoDecEngine(use_wrapper=False, verbose=False)
eng_thrash.config.activescan = 1
eng_thrash.config.activescanrt = -1
eng_thrash.config.activescanorder = 2
eng_thrash.config.mass_diff_c = 1.0033
eng_thrash.predmode = 2

thrash_true_pos = 0
thrash_false_pos = 0
thrash_false_neg = 0

zs_all = {}
zs_nn_only = {}
zs_thrash_only = {}
zs_unmatched = {}

allzs_annotated = {}
allzs_nn = {}
allzs_thrash = {}

ca_seq = 'S[42.010565]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK'
apo_seq = "GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"

ecd_iontypes = ["b", "y", "c", "z"]
etd_iontypes = ['c', 'z']
uvpd_iontypes = ['a','b','c','x','y','z','a+1','x+1','y-1']


#region ToF
tof_spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt")
tof_manual_annotations = read_manual_annotations("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt",
                                             delimiter='\t', data=tof_spectrum)
for p in tof_manual_annotations.peaks:
    if p.z in allzs_annotated:
        allzs_annotated[p.z] += 1
    else:
        allzs_annotated[p.z] = 1

eng_nn.batch_process_spectrum(tof_spectrum, centroided=True)
eng_thrash.batch_process_spectrum(tof_spectrum, centroided=True)

#get_inconsistincies(tof_spectrum, tof_manual_annotations, eng_nn.pks, eng_thrash.pks, maxshift=3, ppmtol=20)

print("IsoDec 8-bin Results on ToF Data...")
shared_nn, unique_annotated_nn, unique_experimental_nn = (
    compare_annotated(tof_manual_annotations, eng_nn.pks, ppmtol=50, maxshift=3))
nn_true_pos, nn_false_pos, nn_false_neg = print_precision_recall(len(shared_nn), len(unique_annotated_nn), len(unique_experimental_nn),
                                                                 nn_true_pos, nn_false_pos, nn_false_neg)
for p in eng_nn.pks.peaks:
    if p.z in allzs_nn:
        allzs_nn[p.z] += 1
    else:
        allzs_nn[p.z] = 1
searchTD(eng_nn.pks, apo_seq, ecd_iontypes, "ppm", 20, verbose=True)

print("IsoDec Thrash Results on ToF Data...")
shared_thrash, unique_annotated_thrash, unique_experimental_thrash = (
    compare_annotated(tof_manual_annotations, eng_thrash.pks, ppmtol=50, maxshift=3))
thrash_true_pos, thrash_false_pos, thrash_false_neg = print_precision_recall(len(shared_thrash), len(unique_annotated_thrash), len(unique_experimental_thrash),
                       thrash_true_pos, thrash_false_pos, thrash_false_neg)

for p in eng_thrash.pks.peaks:
    if p.z in allzs_thrash:
        allzs_thrash[p.z] += 1
    else:
        allzs_thrash[p.z] = 1
searchTD(eng_thrash.pks, apo_seq, ecd_iontypes, "ppm", 20, verbose=True)

#Now get the peaks that are unique to the NN and THRASH
shared, u1, u2 = compare_matchedcollections(shared_nn, shared_thrash, ppmtol=50, maxshift=3, ignorescan=True)
for p in shared.peaks:
    if p.z in zs_all:
        zs_all[p.z] += 1
    else:
        zs_all[p.z] = 1

for p in u1.peaks:
    if p.z in zs_nn_only:
        zs_nn_only[p.z] += 1
    else:
        zs_nn_only[p.z] = 1

for p in u2.peaks:
    if p.z in zs_thrash_only:
        zs_thrash_only[p.z] += 1
    else:
        zs_thrash_only[p.z] = 1

unmatched_annotated, u1, u2 = compare_matchedcollections(unique_annotated_nn, unique_annotated_thrash, ppmtol=50, maxshift=3, ignorescan=True)

for p in unmatched_annotated.peaks:
    if p.z in zs_unmatched:
        zs_unmatched[p.z] += 1
    else:
        zs_unmatched[p.z] = 1
#endregion ToF

eng_nn.pks = MatchedCollection()
eng_thrash.pks = MatchedCollection()

#region ETD
etd_spectrum = np.loadtxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\etd_spectrum.txt")
etd_manual_annotations = read_manual_annotations("Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt",
                                             delimiter=' ')
for p in etd_manual_annotations.peaks:
    if p.z in allzs_annotated:
        allzs_annotated[p.z] += 1
    else:
        allzs_annotated[p.z] = 1

eng_nn.batch_process_spectrum(etd_spectrum, centroided=True)
eng_thrash.batch_process_spectrum(etd_spectrum, centroided=True)

print("IsoDec 8-bin Results on ETD Data...")
shared_nn, unique_annotated_nn, unique_experimental_nn = (
    compare_annotated(etd_manual_annotations, eng_nn.pks, ppmtol=50, maxshift=3))
nn_true_pos, nn_false_pos, nn_false_neg = print_precision_recall(len(shared_nn), len(unique_annotated_nn), len(unique_experimental_nn),
                       nn_true_pos, nn_false_pos, nn_false_neg)
for p in eng_nn.pks.peaks:
    if p.z in allzs_nn:
        allzs_nn[p.z] += 1
    else:
        allzs_nn[p.z] = 1
searchTD(eng_nn.pks, ca_seq, etd_iontypes, "ppm", 20, verbose=True)

print("IsoDec Thrash Results on ETD Data...")
shared_thrash, unique_annotated_thrash, unique_experimental_thrash = (
    compare_annotated(etd_manual_annotations, eng_thrash.pks, ppmtol=50, maxshift=3))
thrash_true_pos, thrash_false_pos, thrash_false_neg = print_precision_recall(len(shared_thrash), len(unique_annotated_thrash), len(unique_experimental_thrash),
                                                                             thrash_true_pos, thrash_false_pos, thrash_false_neg)
for p in eng_thrash.pks.peaks:
    if p.z in allzs_thrash:
        allzs_thrash[p.z] += 1
    else:
        allzs_thrash[p.z] = 1
searchTD(eng_thrash.pks, ca_seq, etd_iontypes, "ppm", 20, verbose=True)

#Now get the peaks that are unique to the NN and THRASH
shared, u1, u2 = compare_matchedcollections(shared_nn, shared_thrash, ppmtol=50, maxshift=3, ignorescan=True)
for p in shared.peaks:
    if p.z in zs_all:
        zs_all[p.z] += 1
    else:
        zs_all[p.z] = 1

for p in u1.peaks:
    if p.z in zs_nn_only:
        zs_nn_only[p.z] += 1
    else:
        zs_nn_only[p.z] = 1

for p in u2.peaks:
    if p.z in zs_thrash_only:
        zs_thrash_only[p.z] += 1
    else:
        zs_thrash_only[p.z] = 1

unmatched_annotated, u1, u2 = compare_matchedcollections(unique_annotated_nn, unique_annotated_thrash, ppmtol=50, maxshift=3, ignorescan=True)

for p in unmatched_annotated.peaks:
    if p.z in zs_unmatched:
        zs_unmatched[p.z] += 1
    else:
        zs_unmatched[p.z] = 1

#endregion ETD

eng_nn.pks = MatchedCollection()
eng_thrash.pks = MatchedCollection()

#region UVPD
uvpd_spectrum = np.loadtxt("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\uvpd_spectrum.txt")
uvpd_manual_annotations = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt",
                                             delimiter='\t')
for p in uvpd_manual_annotations.peaks:
    if p.z in allzs_annotated:
        allzs_annotated[p.z] += 1
    else:
        allzs_annotated[p.z] = 1

eng_nn.batch_process_spectrum(uvpd_spectrum, centroided=True)
eng_thrash.batch_process_spectrum(uvpd_spectrum, centroided=True)

print("IsoDec 8-bin Results on UVPD Data...")
shared_nn, unique_annotated_nn, unique_experimental_nn = (
    compare_annotated(uvpd_manual_annotations, eng_nn.pks, ppmtol=50, maxshift=3))
nn_true_pos, nn_false_pos, nn_false_neg = print_precision_recall(len(shared_nn), len(unique_annotated_nn), len(unique_experimental_nn),
                                                                 nn_true_pos, nn_false_pos, nn_false_neg)
for p in eng_nn.pks.peaks:
    if p.z in allzs_nn:
        allzs_nn[p.z] += 1
    else:
        allzs_nn[p.z] = 1
searchTD(eng_nn.pks, ca_seq, uvpd_iontypes, "ppm", 20, verbose=True)

print("IsoDec Thrash Results on UVPD Data...")
shared_thrash, unique_annotated_thrash, unique_experimental_thrash = (
    compare_annotated(uvpd_manual_annotations, eng_thrash.pks, ppmtol=50, maxshift=3))
thrash_true_pos, thrash_false_pos, thrash_false_neg = print_precision_recall(len(shared_thrash), len(unique_annotated_thrash), len(unique_experimental_thrash),
                       thrash_true_pos, thrash_false_pos, thrash_false_neg)
for p in eng_thrash.pks.peaks:
    if p.z in allzs_thrash:
        allzs_thrash[p.z] += 1
    else:
        allzs_thrash[p.z] = 1
searchTD(eng_thrash.pks, ca_seq, uvpd_iontypes, "ppm", 20, verbose=True)


#Now get the peaks that are unique to the NN and THRASH
shared, u1, u2 = compare_matchedcollections(shared_nn, shared_thrash, ppmtol=50, maxshift=3, ignorescan=True)
for p in shared.peaks:
    if p.z in zs_all:
        zs_all[p.z] += 1
    else:
        zs_all[p.z] = 1

for p in u1.peaks:
    if p.z in zs_nn_only:
        zs_nn_only[p.z] += 1
    else:
        zs_nn_only[p.z] = 1

for p in u2.peaks:
    if p.z in zs_thrash_only:
        zs_thrash_only[p.z] += 1
    else:
        zs_thrash_only[p.z] = 1

unmatched_annotated, u1, u2 = compare_matchedcollections(unique_annotated_nn, unique_annotated_thrash, ppmtol=50, maxshift=3, ignorescan=True)

for p in unmatched_annotated.peaks:
    if p.z in zs_unmatched:
        zs_unmatched[p.z] += 1
    else:
        zs_unmatched[p.z] = 1
#endregion UVPD

for i in range(1,51):
    if i not in allzs_annotated:
        allzs_annotated[i] = 0
    if i not in zs_all:
        zs_all[i] = 0
    if i not in zs_nn_only:
        zs_nn_only[i] = 0
    if i not in zs_thrash_only:
        zs_thrash_only[i] = 0
    if i not in zs_unmatched:
        zs_unmatched[i] = 0


zs_all = {key: zs_all[key] for key in sorted(zs_all)}
zs_nn_only = {key: zs_nn_only[key] for key in sorted(zs_nn_only)}
zs_thrash_only = {key: zs_thrash_only[key] for key in sorted(zs_thrash_only)}
zs_unmatched = {key: zs_unmatched[key] for key in sorted(zs_unmatched)}
allzs_annotated = {key: allzs_annotated[key] for key in sorted(allzs_annotated)}

print("IsoDec 8-bin Results on All Data...")
print_precision_recall(nn_true_pos, nn_false_neg, nn_false_pos, 0,0,0)
print("IsoDec Thrash Results on All Data...")
print_precision_recall(thrash_true_pos, thrash_false_neg, thrash_false_pos, 0,0,0)


if True:
    plt.plot(list(allzs_annotated.keys()), list(allzs_annotated.values()), label="Annotated")
    plt.plot(list(zs_all.keys()), list(zs_all.values()), label="Both Matched")
    plt.plot(list(zs_nn_only.keys()), list(zs_nn_only.values()), label="NN Only")
    plt.plot(list(zs_thrash_only.keys()), list(zs_thrash_only.values()), label="THRASH Only")
    plt.plot(list(zs_unmatched.keys()), list(zs_unmatched.values()), label="Unmatched")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylim(0, 300)
    plt.xticks(np.arange(0, 51, step=5), fontsize=16)
    plt.yticks(np.arange(0, 301, step=100), fontsize=16)
    plt.xlim(1, 50)
    plt.xlabel("Charge", fontsize=24)
    plt.ylabel("Count", fontsize=24)
    plt.legend()
    plt.tight_layout()
    if False:
        plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\zdist_matched_unmatched.png", dpi=300, transparent=True)
        plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\zdist_matched_unmatched.pdf", dpi=300, transparent=True)
    plt.show()


allzs_nn = {key: allzs_nn[key] for key in sorted(allzs_nn)}
allzs_nn_fmatched = {}
for key in allzs_nn:
    total = allzs_nn[key]
    matched = zs_all.get(key, 0) + zs_nn_only.get(key, 0)
    allzs_nn_fmatched[key] = matched / total
allzs_nn_fmatched = {key: allzs_nn_fmatched[key] for key in sorted(allzs_nn_fmatched)}

allzs_thrash = {key: allzs_thrash[key] for key in sorted(allzs_thrash)}
allzs_thrash_fmatched = {}
for key in allzs_thrash:
    total = allzs_thrash[key]
    matched = zs_all.get(key, 0) + zs_thrash_only.get(key, 0)
    allzs_thrash_fmatched[key] = matched / total
allzs_thrash_fmatched = {key: allzs_thrash_fmatched[key] for key in sorted(allzs_thrash_fmatched)}

plt.plot(list(allzs_nn_fmatched.keys()), list(allzs_nn_fmatched.values()), label="NN")
plt.plot(list(allzs_thrash_fmatched.keys()), list(allzs_thrash_fmatched.values()), label="Thrash")
plt.xlim(0, 30)
plt.ylim(0, 1)
plt.xlabel("Charge")
plt.ylabel("Fraction Matched")
plt.legend()
plt.show()