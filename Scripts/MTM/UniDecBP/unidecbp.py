import os
from unidec import engine
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def make_plot(eng, matches, outname=None):
    massdata = eng.data.massdat
    fig = plt.figure()
    plt.subplot(121)
    plt.plot(eng.data.data2[:,0], eng.data.data2[:,1])

    plt.subplot(122)
    plt.plot(massdata[:,0], massdata[:,1])
    plt.xlim(np.amin(eng.pks.masses)-2000, np.amax(eng.pks.masses)+2000)

    for i, match in enumerate(matches):
        p = eng.pks.peaks[i]
        m = p.mass
        h = p.height

        if "Correct" in match:
            color="g"
            plt.text(m, h+5, match)
        elif match == "":
            color = "y"
        else:
            color = "r"
            plt.text(m, h + 5, match)

        plt.plot(m, h, color=color, marker="o")

    plt.show()
    #html_writer.fig_to_html(fig, outname)


tolerance = 50

seqfile = "C:\Data\Wilson_Genentech\sequences.xlsx"
rawfolder = "C:\Data\Wilson_Genentech\Raw Files for Michael Marty"
modfile = "C:\Data\Wilson_Genentech\mods.csv"

os.chdir(rawfolder)

seqdf = pd.read_excel(seqfile)
moddf = pd.read_csv(modfile)

eng = engine.UniDec()

run = False
time_range = [3.9, 4.6]


percs = []
percs2 = []
matches = []

for i, file in enumerate(seqdf["Sample name"]):
    fname = file + ".raw"

    if os.path.isfile(fname):
        print(fname)
        eng.open_file(fname, rawfolder, time_range=time_range)
        eng.config.peakwindow = 15
        eng.config.peakthresh = 0.01
        eng.config.massub = 200000
        if run:
            eng.autorun()
        else:
            eng.unidec_imports(efficiency=True)
            eng.pick_peaks()

        perc, perc2, matcharray, matchstring = UPP_check_peaks(seqdf.iloc[i], eng.pks, tolerance, moddf)
        percs.append(perc)
        percs2.append(perc2)
        matches.append(matchstring)

        make_plot(eng, matcharray, fname[:-4]+".html")
    else:
        print("File not found:", fname, rawfolder)

percs = np.array(percs)
percs2 = np.array(percs2)

seqdf["Correct %"] = percs[:, 0]
seqdf["Incorrect %"] = percs[:, 1]
seqdf["Unmatched %"] = percs[:, 2]
seqdf["Correct % Matched Only"] = percs2[:, 0]
seqdf["Incorrect % Matched Only"] = percs2[:, 0]
seqdf["Matches"] = matches

seqdf.to_excel("results.xlsx")
