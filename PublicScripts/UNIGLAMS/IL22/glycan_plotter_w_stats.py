import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.engine as engine
import unidec.tools as ud
import matplotlib.cm as cm

def make_plot(axs, i, sdata, gdata, hdata, fdata):
    alpha = 0.4
    axs[i].bar(sdata[:, 0], sdata[:, 1], alpha=alpha)
    axs[i].plot(sdata[:, 0], sdata[:, 1], label="Sialic Acid", alpha=alpha)
    axs[i].bar(gdata[:, 0], gdata[:, 1], alpha=alpha)
    axs[i].plot(gdata[:, 0], gdata[:, 1], label="GlcNAc", alpha=alpha)
    axs[i].bar(hdata[:, 0], hdata[:, 1], alpha=alpha)
    axs[i].plot(hdata[:, 0], hdata[:, 1], label="Hexose", alpha=alpha)
    axs[i].bar(fdata[:, 0], fdata[:, 1], alpha=alpha)
    axs[i].plot(fdata[:, 0], fdata[:, 1], label="Fucose", alpha=alpha)
    if i == 0:
        axs[i].legend()
    # axs[i].set_title(pfile[:-4])
    axs[i].set_title(titles[i])
    axs[i].set_xlim(0, 30)
    axs[i].set_ylim(0, 0.8)
    axs[i].set_ylabel("Total Probability")


def make_2D_plot(axs, i, masses, sdata, sdata2d, gdata, gdata2d, hdata, hdata2d, fdata, fdata2d):
    alpha = 0.5
    nlevels = 100
    normalization = cm.colors.Normalize(vmax=0.75, vmin=0.0)

    bluecm = ud.make_alpha_cmap((0,0,1), alpha)
    redcm = ud.make_alpha_cmap((1, 0, 0), alpha)
    greencm = ud.make_alpha_cmap((0, 1, 0), alpha)
    orangecm = ud.make_alpha_cmap((1, 0.65, 0), alpha)

    axs[i].contourf(masses/1000., sdata[:, 0], sdata2d,levels=nlevels, cmap=bluecm, norm=normalization)
    axs[i].contourf(masses / 1000., gdata[:, 0], gdata2d, levels=nlevels, cmap=orangecm, norm=normalization)
    axs[i].contourf(masses / 1000., hdata[:, 0], hdata2d, levels=nlevels, cmap=greencm, norm=normalization)
    axs[i].contourf(masses / 1000., fdata[:, 0], fdata2d, levels=nlevels, cmap=redcm, norm=normalization)
    axs[i].set_xlim(95, 110)
    axs[i].set_ylim(0, 30)
    axs[i].set_ylabel("Total Number")




os.chdir("C:\\Data\\Luis Genentech")

eng = engine.UniDec()

pfile1 = "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_peaks.txt"
pfile2 = "20220113144226FcFusionProtein_SA8_20uScans_peaks.txt"
pfile3 = "211029_Wendy_Isle22FC_4_warmup_PTR2_peaks.txt"

files = [pfile3, pfile2, pfile1]

statsfiles = ["IL22_SA4_tryp_glycan_list_stats.npy", "IL22_SA8_tryp_glycan_list_stats.npy",
              "SA15_Trypsin_glycan_list_stats.npy"]

# statsfiles = ["IL22_SA4_LysC_v2_glycan_list_stats.npy", "IL22_SA8_LysC_glycan_list_stats.npy",
#              "SA15_LysC_glycan_list_stats.npy"]

ofiles = ["SA4_ofile_filtered.dat", "SA8_ofile_filtered.dat",
          "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile_filtered.dat"]

ofiles = ["20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile.dat",
          "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile.dat",
          "20220113162114FcFusionProtein_SA15_20uScans_5500_7000mz_ofile.dat"]

titles = ["SA4", "SA8", "SA15"]

glycomode = True
tolerance = 5
normmode = 1
alts = True

fig, axs = plt.subplots(len(files), figsize=(5, 10))
fig2, axs2 = plt.subplots(len(files), figsize=(5, 10))
fig3, axs3 = plt.subplots(len(files), figsize=(5, 10))
for i, pfile in enumerate(files):
    eng.load_peaks(pfile)

    eng.load_ofile(ofiles[i])

    eng.match(tolerance=tolerance, glyco=glycomode)  # , minsites=6, maxsites=8)
    eng.get_alts(tolerance=tolerance)
    print(eng.matchcounts)

    sindex, hindex, gindex, findex = ud.get_glyco_indexes(eng.olg.oligomerlist)

    statsarray = np.load(statsfiles[i], allow_pickle=True)
    statsarray = [None, None, None, None]

    sdata, sdata2d = eng.get_summed_match_intensities(sindex, normmode=normmode, alts=alts, probarray=statsarray[0],
                                                      get_2D=True)
    gdata, gdata2d = eng.get_summed_match_intensities(gindex, normmode=normmode, alts=alts, probarray=statsarray[1],
                                                      get_2D=True)
    hdata, hdata2d = eng.get_summed_match_intensities(hindex, normmode=normmode, alts=alts, probarray=statsarray[2],
                                                      get_2D=True)
    fdata, fdata2d = eng.get_summed_match_intensities(findex, normmode=normmode, alts=alts, probarray=statsarray[3],
                                                      get_2D=True)
    statsarray = np.load(statsfiles[i], allow_pickle=True)
    make_plot(axs, i, sdata, gdata, hdata, fdata)

    make_plot(axs2, i, statsarray[0], statsarray[1], statsarray[2], statsarray[3])

    make_2D_plot(axs3, i, eng.pks.masses, sdata, sdata2d, gdata, gdata2d, hdata, hdata2d, fdata, fdata2d)
    index = i


axs[index].set_xlabel("Number of Total Units")
fig.tight_layout()
if True:
    outfileheader = pfile[:-4] + "_match_glycomode" + str(glycomode) + "_tol_" + str(tolerance) + "_NativeOnly"
    fig.savefig(outfileheader + ".png")
    fig.savefig(outfileheader + ".pdf")

axs2[index].set_xlabel("Number of Total Units")
fig2.tight_layout()
if True:
    outfileheader = outfileheader + '_bottomUP'
    fig2.savefig(outfileheader + ".png")
    fig2.savefig(outfileheader + ".pdf")

axs3[index].set_xlabel("Mass (kDa)")
fig3.tight_layout()
if True:
    outfileheader = outfileheader + '_2D'
    fig3.savefig(outfileheader + ".png")
    fig3.savefig(outfileheader + ".pdf")

fig.show()
fig2.show()
fig3.show()
