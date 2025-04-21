import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec import engine

'''
eng = unidec.unidec()
eng.open_test_spectrum(masslist=[12000], intlist=[1], noise=0.02, resolution=1000, psfun=2, silent=True, zwidth=0.75)
eng.config.default_decon_params()
eng.config.default_high_res()
eng.config.massbins = 1
eng.config.masslb = 11500
eng.config.massub = 12500
eng.process_data()
eng.run_unidec()
eng.pick_peaks()

exit()
'''

np.random.seed(25)

acolor = "g"
ccolor = "purple"
pcolor = "navy"
abcindex = ord("a")

savepath = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\MSError\\Figure2.pdf"
savepath2 = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\MSError\\Figure2.png"

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
tpad = 10

with plt.xkcd(scale=0.5, length=100, randomness=1):
    plt.figure(figsize=(6.5, 6))
    plt.subplot(311)

    md = []
    for i in range(0, 3):
        eng = engine.UniDec()
        eng.open_test_spectrum(masslist=[12000], intlist=[1], noise=0.03, resolution=1000, psfun=2, silent=True,
                               zwidth=0.75, mz_bin_size=0.1)
        eng.config.default_decon_params()
        eng.config.default_high_res()
        eng.config.massbins = 0.1
        eng.config.masslb = 11500
        eng.config.massub = 12500
        eng.process_data()
        eng.run_unidec()
        eng.pick_peaks()

        mdata = eng.data.massdat
        mdata[:, 1] /= np.amax(mdata[:, 1])
        data = eng.data.rawdata
        data = np.append(data, np.array([[5500, 0]]), axis=0)

        md.append(mdata)

    peaks = ud.peakdetect(data, threshold=0.1)

    pk = eng.pks.peaks[0]
    peaks = pk.mztab[:, :2]
    threshold=0.1
    b1 = peaks[:, 1] > threshold
    peaks = peaks[b1]

    std = np.std(pk.errormasses[b1])
    wstd = ud.weighted_std(pk.errormasses[b1], pk.errorints[b1])
    print(std, wstd)

    plt.plot(data[:, 0], data[:, 1], color=pcolor)

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])

    xoff = [0, 0, 0, 0, -100]
    for i, p in enumerate(peaks):
        ax.annotate('M' + str(i + 1).translate(SUB), xy=(p[0], p[1]), arrowprops=dict(arrowstyle='->'),
                    xytext=(p[0]+xoff[i], p[1] + 0.3), color=ccolor)

    label = chr(abcindex).upper()
    abcindex += 1
    ax.text(np.amax(data[:, 0]), 1.3, label, horizontalalignment="right", verticalalignment="top", color="k",
            fontsize=14)

    mathstring1 = "$Uncertainty = Std. Dev.(M_1, M_2, M_3, M_4, M_5)$ = " + str(np.round(wstd, 2))
    ax.text(np.amax(data[:, 0]) + 0, 0.8, mathstring1, horizontalalignment="right", verticalalignment="top",
            color=ccolor, fontsize=13)

    plt.title("Std. Dev. of Charge States", pad=tpad)
    plt.xlabel("$\it{m/z}$")
    plt.ylabel("Intensity")
    plt.ylim(0, 1.1)
    plt.xlim(np.amin(data[:, 0]), np.amax(data[:, 0]))

    plt.subplot(312)

    plt.plot(mdata[:, 0], mdata[:, 1], color=pcolor)

    fwhm = eng.pks.peaks[0].intervalFWHM

    plt.hlines(0.5, fwhm[0], fwhm[1], linestyles="--", color=acolor)
    # plt.vlines(centroid, 0, 1, linestyles="--", color=ccolor)

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])

    ax.annotate('Uncertainty = FWHM / 2 = ' + str(str(np.round(eng.pks.peaks[0].errorFWHM / 2, 1))), xy=(fwhm[1], 0.5),
                arrowprops=dict(arrowstyle='->'), xytext=(fwhm[1] + 10, 0.7),
                color=acolor, fontsize=13)
    # ax.annotate('Centroid', xy=(centroid, 1), arrowprops=dict(arrowstyle='->'), xytext=(7, 1.2), color=ccolor)
    label = chr(abcindex).upper()
    abcindex += 1
    ax.text(12060, 1.3, label, horizontalalignment="right", verticalalignment="top", color="k",
            fontsize=14)

    plt.title("Peak Width", pad=tpad)
    plt.xlabel("Mass")
    plt.ylabel("Intensity")
    plt.ylim(0, 1.1)
    plt.xlim(11980, 12060)

    plt.subplot(313)

    colors = [pcolor, "blue", "royalblue"]
    apexes = []
    xoff = [-3, 3, 0]
    for i, m in enumerate(md):
        plt.plot(m[:, 0], (m[:, 1] + i * 0.25) / 1.5, color=colors[i])
        apex = ud.localmaxpos(m, 0, 1000000)
        apexes.append(apex)
        maxval = (1 + i * 0.25) / 1.5
        plt.vlines(apex, 0, maxval, linestyles="--", color=colors[i])
        ax = plt.gca()
        ax.annotate('M' + str(i + 1).translate(SUB), xy=(apex, maxval), arrowprops=dict(arrowstyle='->'),
                    xytext=(apex + xoff[i], maxval + 0.3), color=colors[i])

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])

    ax.text(12060, 0.8, '$Uncertainty = Std. Dev.(M_1, M_2, M_3)$ = ' + str(np.round(np.std(apexes), 2)),
            horizontalalignment="right", verticalalignment="top",
            color="blue", fontsize=13)

    label = chr(abcindex).upper()
    abcindex += 1
    ax.text(12060, 1.3, label, horizontalalignment="right", verticalalignment="top", color="k",
            fontsize=14)

    plt.title("Replicates", pad=tpad)
    plt.xlabel("Mass")
    plt.ylabel("Intensity")
    plt.ylim(0, 1.1)
    plt.xlim(11980, 12060)

    plt.tight_layout()
    plt.savefig(savepath)
    plt.savefig(savepath2)
    plt.show()
