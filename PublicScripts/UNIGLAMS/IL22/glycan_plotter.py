import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.engine as engine
import unidec.modules.unidectools as ud
from unidec.modules.isolated_packages.biopolymer_calculator import calc_pep_mass

os.chdir("C:\\Data\\Luis Genentech")

eng = engine.UniDec()

pfile = "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_peaks.txt"
# pfile = "20220113144226FcFusionProtein_SA8_20uScans_peaks.txt"
# pfile = "211029_Wendy_Isle22FC_4_warmup_PTR2_peaks.txt"
# pfile = "20220711_rhEndosialin_CD248_PTCR_peaks.txt"

eng.load_peaks(pfile)

ofile = "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile.dat"
eng.load_ofile(ofile)
print(eng.config.oligomerlist)

statsfile ="SA15_Trypsin_glycan_list_stats.npy"
statsarray = np.load(statsfile, allow_pickle=True)
print(statsarray)

if False:
    seq = "LLRLLLAWAAAGPTLGQDPWAAEPRAACGPSSCYALFPRRRTFLEAWRACRELGGDLATPRTPEEAQRVDSLVGAGPASRLLWIGLQRQARQCQLQRPLRGFTWTTGDQDTAFTNWAQPASGGPCPAQRCVALEASGEHRWLEGSCTLAVDGYLCQFGFEGACPALQDEAGQAGPAVYTTPFHLVSTEFEWLPFGSVAAVQCQAGRGASLLCVKQPEGGVGWSRAGPLCLGTGCSPDNGGCEHECVEEVDGHVSCRCTEGFRLAADGRSCEDPCAQAPCEQQCEPGGPQGYSCHCRLGFRPAEDDPHRCVDTDECQIAGVCQQMCVNYVGGFECYCSEGHELEADGISCSPAGAMGAQASQDLGDELLDDGEDEEDEDEAWKAFNGGWTEMPGILWMEPTQPPDFALAYRPSFPEDREPQIPYPEPTWPPPLSAPRVPYHSSVLSVTRPVVVSATHPTLPSAHQPPVIPATHPALSRDHQIPVIAANYPDLPSAYQPGILSVSHSAQPPAHQPPMISTKYPELFPAHQSPMFPDTRVAGTQTTTHLPGIPPNHAPLVTTLGAQLPPQAPDALVLRTQATQLPIIPTAQPSLTTTSRSPVSPAHQISVPAATQPAALPTLLPSQSPTNQTSPISPTHPHSKAPQIPREDGPSPKLALWLPSPAPTAAPTALGEAGLAEHSQRDDRWLLVALLVPTCVFLVVLLALGIVYCTRCGPHAPNKRITDCYRWVIHAGSKSPTEPMPPRGSLTGVQTCRTSV"
    seqmass = calc_pep_mass(seq)
    print(seqmass)
    eng.config.oligomerlist[5, 1] = str(seqmass)
    eng.config.oligomerlist[5, 2] = "1"
    eng.config.oligomerlist[5, 3] = "1"
    eng.config.oligomerlist[5, 4] = "P"

tolerance = 5

eng.match(tolerance=tolerance, glyco=True)  # , minsites=6, maxsites=8)
eng.get_alts(tolerance=tolerance)
print(eng.matchcounts)
print(eng.altmasses[0], eng.altindexes[0])

sindex, hindex, gindex, findex = ud.get_glyco_indexes(eng.olg.oligomerlist)

normmode = 1
alts = True
statsarray = [None, None, None, None]
sdata = eng.get_summed_match_intensities(sindex, normmode=normmode, alts=alts, probarray=statsarray[0])
gdata = eng.get_summed_match_intensities(gindex, normmode=normmode, alts=alts, probarray=statsarray[1])
hdata = eng.get_summed_match_intensities(hindex, normmode=normmode, alts=alts, probarray=statsarray[2])
fdata = eng.get_summed_match_intensities(findex, normmode=normmode, alts=alts, probarray=statsarray[3])

alpha=0.4
plt.bar(sdata[:, 0], sdata[:, 1], alpha=alpha)
plt.plot(sdata[:, 0], sdata[:, 1], label="Sialic Acid", alpha=alpha)
plt.bar(gdata[:, 0], gdata[:, 1], alpha=alpha)
plt.plot(gdata[:, 0], gdata[:, 1], label="GlcNAc", alpha=alpha)
plt.bar(hdata[:, 0], hdata[:, 1], alpha=alpha)
plt.plot(hdata[:, 0], hdata[:, 1], label="Hexose", alpha=alpha)
plt.bar(fdata[:, 0], fdata[:, 1], alpha=alpha)
plt.plot(fdata[:, 0], fdata[:, 1], label="Fucose", alpha=alpha)

plt.legend()
plt.title(pfile[:-4])
plt.savefig(pfile[:-4] + "_match.png")
plt.show()
