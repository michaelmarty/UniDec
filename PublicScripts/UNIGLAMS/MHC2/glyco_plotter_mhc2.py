import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.engine as engine
import unidec.modules.unidectools as ud
from unidec.modules.isolated_packages.biopolymer_calculator import calc_pep_mass
import matplotlib.cm as cm

def make_2D_plot(axs, masses, sdata, sdata2d, gdata, gdata2d, hdata, hdata2d, fdata, fdata2d):
    alpha = 0.5
    nlevels = 100
    normalization = cm.colors.Normalize(vmax=0.75, vmin=0.0)

    bluecm = ud.make_alpha_cmap((0,0,1), alpha)
    redcm = ud.make_alpha_cmap((1, 0, 0), alpha)
    greencm = ud.make_alpha_cmap((0, 1, 0), alpha)
    orangecm = ud.make_alpha_cmap((1, 0.65, 0), alpha)

    axs.contourf(masses/1000., sdata[:, 0], sdata2d,levels=nlevels, cmap=bluecm, norm=normalization)
    axs.contourf(masses / 1000., gdata[:, 0], gdata2d, levels=nlevels, cmap=orangecm, norm=normalization)
    axs.contourf(masses / 1000., hdata[:, 0], hdata2d, levels=nlevels, cmap=greencm, norm=normalization)
    axs.contourf(masses / 1000., fdata[:, 0], fdata2d, levels=nlevels, cmap=redcm, norm=normalization)
    #axs.set_xlim(95, 110)
    axs.set_ylim(0, 32)
    axs.set_ylabel("Total Number")
    axs.set_xlabel("Mass (kDa)")

os.chdir("C:\\Data\\Luis Genentech\\MHC2")

eng = engine.UniDec()

pfile = "MHC2_DPA_Combined_peaks.txt"
#pfile = "12012022_MHC2DPA_MS2_PTCR_20221201051533_peaks2.txt"

eng.load_peaks(pfile)

ofile = "MHC2_DPA_Combined_ofile.dat"
#ofile = "12012022_MHC2DPA_MS2_PTCR_20221201051533_ofile.dat"
eng.load_ofile(ofile)
print(eng.config.oligomerlist)

statsfile ="MHC2_DPA_tryp_glycan_list_stats.npy"
statsarray = np.load(statsfile, allow_pickle=True)
print(eng.pks.masses)

tolerance = 5

eng.match(tolerance=tolerance, glyco=True)  # , minsites=6, maxsites=8)
eng.get_alts(tolerance=tolerance)
print(eng.matchcounts)

sindex, hindex, gindex, findex = ud.get_glyco_indexes(eng.olg.oligomerlist)

normmode = 1
alts = True
statsarray = [None, None, None, None]
sdata, sdata2d = eng.get_summed_match_intensities(sindex, normmode=normmode, alts=alts, probarray=statsarray[0], get_2D=True)
gdata, gdata2d = eng.get_summed_match_intensities(gindex, normmode=normmode, alts=alts, probarray=statsarray[1], get_2D=True)
hdata, hdata2d = eng.get_summed_match_intensities(hindex, normmode=normmode, alts=alts, probarray=statsarray[2], get_2D=True)
fdata, fdata2d = eng.get_summed_match_intensities(findex, normmode=normmode, alts=alts, probarray=statsarray[3], get_2D=True)

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

fig3 = plt.figure()#figsize=(5, 10))
axs = plt.gca()
make_2D_plot(axs, eng.pks.masses, sdata, sdata2d, gdata, gdata2d, hdata, hdata2d, fdata, fdata2d)
outfileheader = pfile[:-4] + '_2D_nativeonly'
fig3.savefig(outfileheader + ".png")
fig3.savefig(outfileheader + ".pdf")

plt.show()
