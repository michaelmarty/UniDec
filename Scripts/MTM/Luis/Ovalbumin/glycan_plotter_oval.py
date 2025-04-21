import matplotlib.pyplot as plt
import os
from unidec import engine
import unidec.tools as ud

os.chdir("C:\\Data\\Luis Genentech\\Ovalbumin")

eng = engine.UniDec()

pfile = "20220805_ovalbumin_MS2_3600-3900_PTR20_Iso20_step10_2_peaks.txt"

eng.load_peaks(pfile)

ofile = "20220805_ovalbumin_MS2_3600-3900_PTR20_Iso20_step10_2_ofile.dat"
eng.load_ofile(ofile)
print(eng.config.oligomerlist)

'''
statsfile ="MHC2_DPA_tryp_glycan_list_stats.npy"
statsarray = np.load(statsfile, allow_pickle=True)'''
print(eng.pks.masses)

tolerance = 5

eng.match(tolerance=tolerance, glyco=True)  # , minsites=6, maxsites=8)
eng.get_alts(tolerance=tolerance)
print(eng.matchcounts)

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
if statsarray[0] is None:
    plt.savefig(pfile[:-4] + "_match_native_only.png")
    plt.savefig(pfile[:-4] + "_match_native_only.pdf")
else:
    plt.savefig(pfile[:-4] + "_match_stats.png")
    plt.savefig(pfile[:-4] + "_match_stats.pdf")
plt.show()
