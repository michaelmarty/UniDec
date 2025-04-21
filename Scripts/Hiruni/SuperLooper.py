import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import matplotlib
import scipy.stats as stats
import statsmodels.api as sm

matplotlib.use('WxAgg')
# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 1


def get_mass_data(newpath):
    files = []
    paths = []
    datalist = []
    for path in glob.glob(f'{newpath}/*/**/', recursive=True):
        mass_file = [f for f in os.listdir(path) if "_mass.txt" in f]
        if len(mass_file) > 0:
            mass_file = os.path.join(path, mass_file[0])
            files.append(mass_file)
            paths.append(path)
            data = np.loadtxt(mass_file)
            data[:, 1] = data[:, 1] / np.amax(data[:, 1])
            data[:, 0] = data[:, 0] / 1000.
            datalist.append(data)
    if len(files) < 3:
        print("Failed at", newpath)
    return datalist


def get_int_data(topdir, fname):
    files = []
    for path in glob.glob(f'{topdir}/*/**/', recursive=True):
        testpath = os.path.join(path, fname)
        if os.path.isfile(testpath):
            files.append(testpath)
    # print(files)
    if len(files) < 3:
        print("Failed at t=", t)
    intdat = []
    for f in files:
        data = np.loadtxt(f)
        intdat.append(data)
    return np.array(intdat)


def get_kd(mutant, wt, n=0, l=1, T=37):
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    TK = T + 273.15
    deltag = -R * TK * np.log(kd)
    return kd, deltag, deltag / l, n, l, n + l


def get_kdg(data, temp, maxl=7):
    bm = data[:, 1] == 0
    mutant = data[bm][:, 2]
    wt = data[np.logical_not(bm)][:, 2]

    wt = wt[:maxl]
    mutant = mutant[:maxl]

    wt = wt / np.amax(data[:, 2])
    mutant = mutant / np.amax(data[:, 2])

    kdlist = []
    for n in range(0, maxl - 1):
        value = get_kd(mutant, wt, n, 1, T=temp)
        kdlist.append(value)
    return wt, mutant, kdlist


def get_avg_kdg(datalist, temp, maxl=7):
    wts = []
    mts = []
    kds = []
    for i in range(len(datalist)):
        data = datalist[i]
        wt, mutant, kdlist = get_kdg(data, temp, maxl)
        kds.append(kdlist)
        wts.append(wt)
        mts.append(mutant)
    kds = np.array(kds)
    wts = np.array(wts)
    mts = np.array(mts)

    wavg = np.mean(wts, axis=0)
    mavg = np.mean(mts, axis=0)

    wstd = np.std(wts, axis=0)
    mstd = np.std(mts, axis=0)

    kavg = np.mean(kds[:, :, 0], axis=0)
    kstd = np.std(kds[:, :, 0], axis=0)

    dgs = []
    for n in range(0, maxl - 1):
        avgdg = np.mean(kds[:, n, 2])
        sddg = np.std(kds[:, n, 2])
        dgvals = list(kds[:, n, 2])
        dgvals.append(avgdg)
        dgvals.append(sddg)
        dgs.append(dgvals)
    return wavg, wstd, mavg, mstd, np.array(dgs), wts, mts, kavg, kstd


def make_mass_plot(ax, massdatalist, p, maxl, wavg, wstd, mavg, mstd, t, e1, e2, edges):
    wmass = (p[0] + p[2] + p[1] * np.arange(0, maxl)) / 1000.
    mmass = (p[0] + p[1] * np.arange(0, maxl)) / 1000.
    for i in range(len(massdatalist)):
        masses = massdatalist[i]
        ax.plot(masses[:, 0], masses[:, 1], color="dimgrey")
        # ax.plot(wmass, wts[i], marker="o", linestyle=" ", color="b")
        # ax.plot(mmass, mts[i], marker="o", linestyle=" ", color="g")

    # ax.bar(wmass, wavg, width=0.2, yerr=wstd, color="k", capsize=5)
    # ax.bar(mmass, mavg, width=0.2, yerr=mstd, color="k", capsize=5)
    ax.errorbar(wmass, wavg, marker="o", linestyle="--", yerr=wstd, color="b", capsize=5)
    ax.errorbar(mmass, mavg, marker="o", linestyle="--", yerr=mstd, color="g", capsize=5)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(0, 1)
    if ed == 0:
        ax.set_title(str(t) + " °C")
    if ed == len(edges) - 1:
        ax.set_xlabel("Mass (kDa)")
    else:
        ax.xaxis.set_ticklabels([])
    if tind == 0:
        ax.set_ylabel(e1 + " " + e2)
        ax.set_yticks([0, 0.5, 1], labels=["0", "%", "100"])
    else:
        ax.set_yticks([0, 0.5, 1], labels=["", "", ""])


def vant_hoff(avgstd, alldata, maxl):
    for l in range(maxl - 1):
        b1 = alldata[:, 0] == l
        data = alldata[b1]

        x = 1 / (data[:, 1] + 273.15)
        y = -data[:, 2] * x / R
        model = sm.OLS(y, sm.add_constant(x))
        results = model.fit()

        df = len(x) - 2
        tval = abs(stats.t.ppf(0.05 / 2, df))
        mH = -results.params[1] * R
        mS = results.params[0] * R
        sH = results.bse[1] * R * tval
        sS = results.bse[0] * R * tval

        preds = results.get_prediction().summary_frame(0.05)
        diffs = (preds["mean_ci_upper"].to_numpy() - preds["mean_ci_lower"].to_numpy())[6]
        dGAvg = (preds["mean"][6] * -R) / x[6]
        dGStd = (diffs / 2. * -R) / x[6]

        mTS = -298.15 * mS
        sTS = -298.15 * sS
        print('L=', l, "H=", mH, "+/-", sH, "S=", mS, "+/-", sS, "TS=", mTS, "+/-", sTS)

        b1 = avgstd[:, 0] == l
        avgstd[b1, 4] = mH
        avgstd[b1, 5] = mS
        avgstd[b1, 6] = mTS
        avgstd[b1, 7] = sH
        avgstd[b1, 8] = sS
        avgstd[b1, 9] = sTS
        avgstd[b1, 10] = dGAvg
        avgstd[b1, 11] = dGStd
    return avgstd


def make_k_plot(ax, kavg, kstd, maxl, edges, t, ed, tind, basis=1, yrange=[0.6, 1.2], cnorm=0.3):
    x = np.arange(1, maxl)

    tval = abs(stats.t.ppf(0.05 / 2, 2))
    failcolor = "grey"
    if basis == 1:
        cmap = plt.cm.bwr_r
    else:
        cmap = plt.cm.bwr
    norm = matplotlib.colors.TwoSlopeNorm(basis, vmin=basis - cnorm, vmax=basis + cnorm)
    colors = []
    for i in x - 1:
        diff = np.abs(basis - kavg[i])
        ci = kstd[i] * tval / np.sqrt(3)
        if diff > ci:
            color = cmap(norm(kavg[i]))
            colors.append(color)
        else:
            colors.append(failcolor)

    ax.bar(x, kavg, yerr=kstd, color=colors)
    ax.plot(x, x * 0 + basis, linestyle="--", color="grey")
    ax.set_ylim(yrange[0], yrange[1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if ed == 0:
        ax.set_title(str(t) + " °C")
    if ed == len(edges) - 1:
        ax.set_xlabel("# Lipids")
        ax.set_xticks([1, 2, 3, 4, 5, 6])  # , labels=["0", "%", "100"])
    else:
        ax.xaxis.set_ticklabels([])
    if tind == 0:
        ax.set_ylabel(e1 + " " + e2)


def make_subplot(ax, x, y, yerr, basis=0, cnorm=0.75, yrange=[-5, 5]):
    failcolor = "grey"
    cmap = plt.cm.bwr
    norm = matplotlib.colors.TwoSlopeNorm(0, vmin=basis - cnorm, vmax=basis + cnorm)
    colors = []
    for i in x - 1:
        diff = np.abs(basis - y[i])
        if diff > np.abs(yerr[i]):
            color = cmap(norm(y[i]))
            colors.append(color)
        else:
            colors.append(failcolor)
    ax.bar(x, y, yerr=yerr, color=colors)
    ax.plot(x, x * 0 + basis, linestyle="--", color="grey")
    ax.set_ylim(yrange[0], yrange[1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def make_tsh_plot(ax, avgstd, maxl, edges, ed, e1, e2, cnorm=0.75, basis=0):
    data = avgstd[:maxl - 1]
    mH = data[:, 4]
    #mS = data[:, 5]
    mTS = data[:, 6]
    sH = data[:, 7]
    #sS = data[:, 8]
    sTS = data[:, 9]
    mDG = data[:, 10]
    sDG = data[:, 11]

    x = np.arange(1, maxl)
    make_subplot(ax[0], x, mDG, sDG, yrange=[-0.5, 1.5])
    make_subplot(ax[1], x, mH, sH)
    make_subplot(ax[2], x, mTS, sTS)
    if ed == 0:
        ax[0].set_title("ΔΔG")
        ax[1].set_title("ΔΔH")
        ax[2].set_title("-TΔΔS")

    for a in ax:
        if ed == len(edges) - 1:
            a.set_xlabel("# Lipids")
            a.set_xticks([1, 2, 3, 4, 5, 6])  # , labels=["0", "%", "100"])
        else:
            a.xaxis.set_ticklabels([])
    ax[0].set_ylabel(e1 + " " + e2)


R = 8.31446261815324 / 1000.
units = "kJ/mol"

fname = "Extract_grid_2D_Extract2.txt"
globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final\\"

os.chdir(globaldir)
edges = [("WT", "R224A"), ("WT", "R75A"), ("WT", "DM"), ("R224A", "DM"), ("R75A", "DM")]
temps = [15, 20, 25, 30, 35, 40]
if False:
    edges = edges[:2]
    # temps = temps[:2]
width = 100
method = 6
params1 = [98590, 1400, 350, 0, 7, 0, 1, width, method]
params2 = [98260, 1400, 680, 0, 7, 0, 1, width, method]
params3 = [98260, 1400, 300, 0, 7, 0, 1, width, method]
params = [params1, params1, params2, params3, params3]
maxl = 5

fig, axs = plt.subplots(ncols=len(temps), nrows=len(edges), figsize=(12, 9),
                        constrained_layout=True)
fig2, axs2 = plt.subplots(ncols=len(temps), nrows=len(edges), figsize=(12, 9),
                          constrained_layout=True)
fig3, axs3 = plt.subplots(ncols=len(temps), nrows=len(edges), figsize=(12, 9),
                          constrained_layout=True)
fig4, axs4 = plt.subplots(ncols=3, nrows=len(edges), figsize=(12, 9),
                          constrained_layout=True)
for ed in range(len(edges)):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    pathtop = os.path.join(globaldir, dirname)
    print(os.path.isdir(pathtop), pathtop)
    out = True
    p = params[ed]

    alldata = []
    avgstd = []
    for tind, t in enumerate(temps):
        newpath = os.path.join(pathtop, e1 + "_" + e2 + "_" + str(t) + "C")

        massdatalist = get_mass_data(newpath)
        intdatalist = get_int_data(newpath, fname)
        wavg, wstd, mavg, mstd, dgs, wts, mts, kavg, kstd = get_avg_kdg(intdatalist, t, maxl=maxl)

        gavg = dgs[:, 3]
        gstd = dgs[:, 4]

        ax = axs[ed, tind]
        ax2 = axs2[ed, tind]
        ax3 = axs3[ed, tind]
        make_mass_plot(ax, massdatalist, p, maxl, wavg, wstd, mavg, mstd, t, e1, e2, edges)
        make_k_plot(ax2, kavg, kstd, maxl, edges, t, ed,tind)
        make_k_plot(ax3, gavg, gstd, maxl, edges, t, ed, tind, yrange=[-0.5, 1.5], basis=0, cnorm=0.75)
        for i, d in enumerate(dgs):
            avgstd.append([i, t, d[-2], d[-1], 0, 0, 0, 0, 0, 0, 0, 0])
            for j, x in enumerate(d[:-2]):
                alldata.append([i, t, x, j])
    alldata = np.array(alldata)
    avgstd = np.array(avgstd)
    vant_hoff(avgstd, alldata, maxl)

    ax4 = axs4[ed]
    make_tsh_plot(ax4, avgstd, maxl, edges, ed, e1, e2)

plt.show()
