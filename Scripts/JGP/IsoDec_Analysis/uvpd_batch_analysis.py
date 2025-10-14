from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.match import *
from Scripts.MTM.proteomicstools import *
import os
import time
from brokenaxes import brokenaxes
from unidec.UniDecImporter.ImporterFactory import ImporterFactory


os.chdir("Z:\\Group Share\\JGP\\JPST001885")

ca_seq = 'S[42.010565]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK'
ion_types = ['a','b','c','x','y','z','a+1','x+1','y-1']

eng = IsoDecEngine(phaseres=8, use_wrapper=True)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

coverages_8bit = []
percs_matched_8bit = []

all_files = ["10uMCA25+UVPDLPT1,5mJ1pFullP.raw", "10uMCA25+UVPDLPT1,5mJ1pFullP_2.raw", "10uMCA25+UVPDLPT1,5mJ1pFullP_3.raw"]

timestart = time.perf_counter()
for file in all_files:
    print("Processing", file)
    reader = ImporterFactory.create_importer(file)
    spectrum = reader.grab_centroid_data(1)
    eng.batch_process_spectrum(spectrum, centroided=True, window=80)
    cov, percent_matched = searchTD(eng.pks.peaks, ca_seq, ion_types, "ppm", 20.0)
    coverages_8bit.append(cov)
    percs_matched_8bit.append(percent_matched)
    eng.pks = MatchedCollection()

print("IsoDec-8bit Time to process:", time.perf_counter()-timestart)
print("IsoDec-8bit Average Coverage:", np.mean(coverages_8bit))
print("IsoDec-8bit Standard Deviation of Coverage:", np.std(coverages_8bit))
print("IsoDec-8bit Average Percent Matched:", np.mean(percs_matched_8bit))
print("IsoDec-8bit Standard Deviation of Percent Matched:", np.std(percs_matched_8bit))

eng = IsoDecEngine(phaseres=4, use_wrapper=True)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

coverages_4bit = []
percs_matched_4bit = []

timestart = time.perf_counter()
for file in all_files:
    print("Processing", file)
    reader = ImporterFactory.create_importer(file)
    spectrum = reader.grab_centroid_data(1)
    eng.batch_process_spectrum(spectrum, centroided=True, window=80)
    cov, percent_matched = searchTD(eng.pks.peaks, ca_seq, ion_types, "ppm", 20.0)
    coverages_4bit.append(cov)
    percs_matched_4bit.append(percent_matched)
    eng.pks = MatchedCollection()

print("IsoDec-4bit Time to process:", time.perf_counter()-timestart)
print("IsoDec-4bit Average Coverage:", np.mean(coverages_4bit))
print("IsoDec-4bit Standard Deviation of Coverage:", np.std(coverages_4bit))
print("IsoDec-4bit Average Percent Matched:", np.mean(percs_matched_4bit))
print("IsoDec-4bit Standard Deviation of Percent Matched:", np.std(percs_matched_4bit))

topfd_msalign_files = ["Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_20241001_094030_TopFD_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP_ms2.msalign",
                       "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_2_20241001_094112_TopFD_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP_2_ms2.msalign",
                       "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_3_20241001_094120_TopFD_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP_3_ms2.msalign"]

top_fd_coverages = []
top_fd_percs_matched = []

for file in topfd_msalign_files:
    uvpd_mc = read_msalign_to_matchedcollection(file, data=None)
    cov, percent_matched = searchTD(uvpd_mc, ca_seq, ion_types, "ppm", 20.0)
    top_fd_coverages.append(cov)
    top_fd_percs_matched.append(percent_matched)

print("TopFD Average Coverage:", np.mean(top_fd_coverages))
print("TopFD Standard Deviation of Coverage:", np.std(top_fd_coverages))
print("TopFD Average Percent Matched:", np.mean(top_fd_percs_matched))
print("TopFD Standard Deviation of Percent Matched:", np.std(top_fd_percs_matched))

ethrash_msalign_files = ["Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_20240930_165240_eTHRASH_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP.msalign",
                         "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_2_20240930_170120_eTHRASH_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP_2.msalign",
                         "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\MASH_Workflows\\10uMCA25+UVPDLPT1,5mJ1pFullP_3_20240930_171136_eTHRASH_1_1\\10uMCA25+UVPDLPT1,5mJ1pFullP_3.msalign"]

ethrash_coverages = []
ethrash_percs_matched = []

for file in ethrash_msalign_files:
    uvpd_mc = read_msalign_to_matchedcollection(file, data=None)
    cov, percent_matched = searchTD(uvpd_mc, ca_seq, ion_types, "ppm", 20.0)
    ethrash_coverages.append(cov)
    ethrash_percs_matched.append(percent_matched)

print("eTHRASH Average Coverage:", np.mean(ethrash_coverages))
print("eTHRASH Standard Deviation of Coverage:", np.std(ethrash_coverages))
print("eTHRASH Average Percent Matched:", np.mean(ethrash_percs_matched))
print("eTHRASH Standard Deviation of Percent Matched:", np.std(ethrash_percs_matched))


flashdeconv_msalign_files = ["Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\FLASHDeconv_msalign\\10uMCA25+UVPDLPT1,5mJ1pFullP_ms2.msalign",
                             "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\FLASHDeconv_msalign\\10uMCA25+UVPDLPT1,5mJ1pFullP_2_ms2.msalign",
                             "Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\FLASHDeconv_msalign\\10uMCA25+UVPDLPT1,5mJ1pFullP_3_ms2.msalign"]

flashdeconv_coverages = []
flashdeconv_percs_matched = []

for file in flashdeconv_msalign_files:
    uvpd_mc = read_msalign_to_matchedcollection(file, data=None)
    cov, percent_matched = searchTD(uvpd_mc, ca_seq, ion_types, "ppm", 20.0)
    flashdeconv_coverages.append(cov)
    flashdeconv_percs_matched.append(percent_matched)

print("FLASHDeconv Average Coverage:", np.mean(flashdeconv_coverages))
print("FLASHDeconv Standard Deviation of Coverage:", np.std(flashdeconv_coverages))
print("FLASHDeconv Average Percent Matched:", np.mean(flashdeconv_percs_matched))
print("FLASHDeconv Standard Deviation of Percent Matched:", np.std(flashdeconv_percs_matched))

#Plot the times to process
decon_engines = ["ID-4", "ID-8", "FLD", "TFD", "eT"]
colors = ["r", "darkred", "b", "g", "orange"]
times = [12.43,12.5,19.61,43,1400]
fig = plt.figure()
bax = brokenaxes(ylims=((0,70),(1345,1415)), hspace=.05)
bax.bar(decon_engines, times, color=colors, edgecolor='none')
bax.set_ylabel("Time (seconds)", labelpad=40, fontsize=20)
bax.tick_params(axis='x', labelsize=18)
bax.tick_params(axis='y', labelsize=14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
bax.standardize_ticks(ybase=25)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_decon_times.pdf", dpi=300)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_decon_times.png", dpi=300)
plt.show()

#Plot the coverages
mean_coverages = [np.mean(coverages_4bit), np.mean(coverages_8bit), np.mean(flashdeconv_coverages), np.mean(top_fd_coverages), np.mean(ethrash_coverages)]
coverage_stddevs = [np.std(coverages_4bit), np.std(coverages_8bit), np.std(flashdeconv_coverages), np.std(top_fd_coverages), np.std(ethrash_coverages)]
plt.bar(decon_engines, mean_coverages, yerr=coverage_stddevs, color=colors, capsize=3, edgecolor='none')
plt.ylabel("Sequence Coverage (%)", fontsize=20)
plt.ylim(0,100)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.yticks(np.arange(0, 101, 25), fontsize=14)
plt.xticks(fontsize=18)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_coverages.pdf", dpi=300)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_coverages.png", dpi=300)
plt.show()

#Plot the percent matched
mean_percs_matched = [np.mean(percs_matched_4bit), np.mean(percs_matched_8bit), np.mean(flashdeconv_percs_matched),
                      np.mean(top_fd_percs_matched), np.mean(ethrash_percs_matched)]

perc_stddevs = [np.std(percs_matched_4bit), np.std(percs_matched_8bit), np.std(flashdeconv_percs_matched),
                np.std(top_fd_percs_matched), np.std(ethrash_percs_matched)]

plt.bar(decon_engines, mean_percs_matched, yerr=perc_stddevs, color=colors, capsize=3, edgecolor='none')
plt.ylabel("Fragments Matched (%)", fontsize=20)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.yticks(np.arange(0, 101, 25), fontsize=14)
plt.xticks(fontsize=18)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_pmatched.pdf", dpi=300)
plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\uvpd_pmatched.png", dpi=300)
plt.show()




