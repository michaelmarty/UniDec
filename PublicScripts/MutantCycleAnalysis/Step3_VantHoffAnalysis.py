# Script to perform the Van't Hoff analysis on the extracted peak areas
# DDG values are calculated with confidence intervals based on the fit.

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.75
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 16
rcParams['font.sans-serif'] = "Arial"

# Some parameters and equations
R = 8.31446261815324 / 1000.
units = "kJ/mol"


def Kfunc(t, h, s):
    return np.exp((-(h / (R * t)) + s / R))

#Selecting the file
topdir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final2\\"

os.chdir(topdir)
# Edges to consider
edges = [("WT", "W14A")]
for ed in range(0, 1):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    # if e2 is "DM":
    # e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    path = os.path.join(topdir, dirname)
    print(os.path.isdir(path), path)
    out = True

    #Setting the Temperature
    temps = [15, 20, 25, 30, 35, 40]
    maxn = 6
    temp_val = 25

    alldata = []
    avgstd = []
    # Read in the data
    for t in temps:
        newpath = os.path.join(path, e1 + "_" + e2 + "_" + str(t) + "C")
        os.chdir(newpath)
        data = np.loadtxt("DG_resultsN2_Fixed.txt")

        for i, d in enumerate(data):
            avgstd.append([i, t, d[-2], d[-1], 0, 0, 0, 0, 0, 0, 0, 0])
            for j, x in enumerate(d[:-2]):
                alldata.append([i, t, x, j])
    os.chdir(path)
    alldata = np.array(alldata)
    avgstd = np.array(avgstd)

    # Perform the fits
    fig = plt.figure()
    fig.suptitle(dirname)
    for l in range(maxn):
        plt.subplot(round(maxn / 2.) + 1, 2, l + 1)
        lipidn = l
        b1 = alldata[:, 0] == lipidn
        data = alldata[b1]
        b1 = avgstd[:, 0] == l
        meanstd = avgstd[b1]

        # Put data in the right form
        x = 1 / (data[:, 1] + 273.15)
        y = -data[:, 2] * x / R
        df = len(x) - 2
        tval = abs(stats.t.ppf(0.05 / 2, df))

        # Actual fit here
        lr = stats.linregress(x, y)
        fit = [lr.slope, lr.intercept]
        errors = np.array([lr.stderr, lr.intercept_stderr])
        print(fit, errors, errors * tval)

        # Make plots
        plt.scatter(x, y)
        plt.plot(x, x * fit[0] + fit[1])
        x2 = 1 / (meanstd[:, 1] + 273.15)
        y2 = -meanstd[:, 2] * x2 / R
        e2 = -meanstd[:, 3] * x2 / R
        plt.errorbar(x2, y2, np.abs(e2))
        plt.xlabel("1/T (K)")
        plt.ylabel("ln K")
        # plt.ylim(-0.5, 0.5)
        plt.title(str(l))

        # Extract the DDH and DDS values
        mH = -fit[0] * R
        mS = fit[1] * R
        sH = errors[0] * R * tval
        sS = errors[1] * R * tval

        # Get confidence intervals
        model = sm.OLS(y, sm.add_constant(x))
        results = model.fit()
        # print(results.params[::-1]* R, results.bse[::-1]* R*tval)
        #print(results.rsquared_adj)
        #print(results.conf_int(0.05)* R)
        #print(results.pvalues)
        #print(results.t_test([0, 1]))
        #print(results.t_test([1, 0]))

        # Find the closest temperature to the desired temperature
        tk = temp_val + 273.15
        xk = 1 / tk
        closeindex = np.argmin(np.abs(x - xk))
        t_number = closeindex

        preds = results.get_prediction().summary_frame(0.05)
        diffs = (preds["mean_ci_upper"].to_numpy() - preds["mean_ci_lower"].to_numpy())[t_number]
        dGAvg = (preds["mean"][t_number] * -R) / x[t_number]
        dGStd = (diffs / 2. * -R) / x[t_number]
        # print(dGAvg, dGStd)
        # plt.scatter(x, y)
        # plt.plot(x, preds["mean"])
        # plt.plot(x, preds["mean_ci_upper"])
        # plt.plot(x, preds["mean_ci_lower"])
        # plt.show()
        # exit()

        mTS = -tk * mS
        sTS = -tk * sS
        # dGAvg = np.mean(data[:,2])
        # dGStd = np.std(data[:,2])
        print('L=', l, "H=", mH, "+/-", sH, "S=", mS, "+/-", sS, "TS=", mTS, "+/-", sTS)

        b1 = avgstd[:, 0] == lipidn
        avgstd[b1, 4] = mH
        avgstd[b1, 5] = mS
        avgstd[b1, 6] = mTS
        avgstd[b1, 7] = sH
        avgstd[b1, 8] = sS
        avgstd[b1, 9] = sTS
        avgstd[b1, 10] = dGAvg
        avgstd[b1, 11] = dGStd

    # Save Results
    if out:
        np.savetxt("AnalysisResults2.csv", avgstd, delimiter=",", header="Lnum,Temp,G,GE,H,S,TS,HE,SE,TSE,GA,GAE")
    # plt.tight_layout()
    # plt.show()

    # Make more plots
    plt.figure()

    for l in range(maxn):
        index = 0
        plt.subplot(round(maxn / 2.) + 1, 2, l + 1)
        b1 = avgstd[:, 0] == l
        meanstd = avgstd[b1]

        ts = meanstd[:, 1]
        for i, t in enumerate(ts):
            c = cm.spring((t - np.amin(ts)) / np.amax(ts))
            plt.bar(index, meanstd[i, 2], yerr=meanstd[i, 3], color=c, tick_label=str(t))
            index += 1
        plt.bar(index, meanstd[i, 4], yerr=meanstd[i, 7], color="g", tick_label="dH")
        index += 1
        plt.bar(index, meanstd[i, 6], yerr=np.abs(meanstd[i, 9]), color="b", tick_label="dTS")
        index += 1

        plt.bar(index, meanstd[i, 10], yerr=np.abs(meanstd[i, 11]), color="y", tick_label="avgdG")
        index += 1
        print("Avg dG=", meanstd[i, 10], "+/-", meanstd[i, 11])
        plt.title("Lipid " + str(l + 1))

        plt.gca().set_ylim(-6, 6)

    if out:
        np.savetxt("AnalysisResults.csv", meanstd, delimiter=",")
        plt.savefig(os.path.join(topdir, "VH_Plot1.pdf"))
    # plt.show()

    plt.figure()
    x = []
    y = []
    e = []
    for l in range(maxn):
        b1 = avgstd[:, 0] == l
        b2 = avgstd[:, 1] == temp_val
        b3 = np.logical_and(b1, b2)
        meanstd = avgstd[b3][0]
        x.append(l + 1)
        y.append(meanstd[10])  # 10 # 2 means avg ddG from data and 10 means avg ddG from fit
        e.append(meanstd[11])  # 11 # or 3 for std ddG from data and 11 for std ddG from fit

    plt.bar(x, y, yerr=np.abs(e))  # , tick_labels=x+1)
    plt.gca().set_ylim(-0.5, 1.0)
    plt.yticks([-0.5, 0, 0.5, 1.0])
    plt.xlabel("Number of Lipids")
    plt.ylabel("Delta G")
    if out:
        plt.savefig(os.path.join(topdir, "VH_Plot2.pdf"))
    # plt.show()

    x3 = []
    y3 = []
    plt.figure()
    for i, t in enumerate(temps):
        # plt.subplot(2, round(len(temps) / 2.), i + 1)

        b1 = avgstd[:, 1] == t
        meanstd = avgstd[b1]

        x2 = meanstd[:, 0]
        y2 = meanstd[:, 2]
        e2 = meanstd[:, 3]
        fit = np.polyfit(x2, y2, 1)
        # plt.plot(x2, x2 * fit[0] + fit[1])
        plt.plot(x2, x2 * 0)

        x3.append(t)
        y3.append(fit[0])
        plt.errorbar(x2, y2, e2, label=str(t))
        plt.xlabel("Number of Lipids")
        plt.ylabel("Delta G")
        # plt.ylim(-1.25, 0.5)
        # plt.title(str(t))
    # plt.tight_layout()
    plt.legend()
    if out:
        plt.savefig(os.path.join(topdir, "VH_Plot3.pdf"))
plt.show()
