"""
Copyright 2022 Michael Marty, University of Arizona

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import subprocess
from copy import deepcopy


def get_fit(x, y, remove_outliers=True, tolerance_multiple=2):
    x = deepcopy(np.array(x))
    y = deepcopy(np.array(y))
    sortindexes = x.argsort()
    x = x[sortindexes]
    y = y[sortindexes]

    fit = np.polyfit(x, y, 1)

    if remove_outliers:
        fitdat = fit[0] * x + fit[1]
        diffs = np.abs(y - fitdat)
        meansd = np.median(diffs)
        b1 = diffs < meansd * tolerance_multiple

        x = x[b1]
        y = y[b1]
        fit = np.polyfit(x, y, 1)
        fitdat = fit[0] * x + fit[1]
    else:
        fitdat = fit[0] * x + fit[1]

    return fit, np.transpose([x, fitdat]), x, y


def write_sonar(calstring, path="_sonar.inf"):
    string = "[VERSION#1]\n\n" \
             "[DATE]\n" \
             "Acquired Date=11-Jul-2022\n" \
             "Acquired Time=05:56:08\n" \
             "[INSTRUMENT]\n" \
             "Quad Resolution=15.0\n" \
             "Collision Cell Pressure=3.0\n" \
             "Collision Cell Gradient=4.0\n" \
             "Ion Energy=0.2\n" \
             "Post Filter=2\n" \
             "Resolution Mode=Resolution\n" \
             "Pusher Period=85\n" \
             "[SONAR#1]\n" \
             "Scan Time=0.5\n" \
             "Start M/Z=150.0\n" \
             "End M/Z=2000.0\n" \
             "Quadrupole Window=20.0\n"
    text_file = open(path, "wt")
    n = text_file.write(string + calstring)
    text_file.close()


def gen_cal(mz, fit, fit1, fit2):
    y = fit[0] * mz + fit[1]
    y1 = fit1[0] * mz + fit1[1]
    y2 = fit2[0] * mz + fit2[1]
    return "Calibration=" + str(mz) + "," + str(y1) + "," + str(y) + "," + str(y2) + "\n"


apexpath = "C:\\MassLynx\\Apex3D\\Apex3D64.exe"

dpath = "Z:\Group Share\\Melanie Odenkirk\\220711_Ecoli_lipid_LCstepped_SONAR.raw"
dpath = "Z:\Group Share\\Melanie Odenkirk\\220810_Ecoli_ND_collection_methods_tes\\220810_EcoliLipid_d7_SONAR_Res_Neg_run1.raw"
dpath = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\\220808_Ecoli_Lightsplash_SONAR_Res_Neg.raw"
os.chdir(dpath)

header = os.path.split(dpath)[1]
header = os.path.splitext(header)[0]
csvpath = os.path.join(dpath, header + "_Apex3DIons.csv")
print(csvpath)


if True:
    call = [apexpath, "-pRawDirName", dpath, "-outputDirName", dpath, "-outputUserDirName", dpath,
            "-lockMassZ1", str(554.262), "-lockmassToleranceAMU", str(0.25), "-writeBinary", "0", "-bCSVOutput", "1",
            "-leThresholdCounts", str(100), "-startingRTMin", str(0), "-endingRTMin", str(20), "-function", str(1)]

    out = subprocess.run(call)

df = pd.read_csv(csvpath)
b1 = df["m_z"] > 400
df = df[b1]

xinit = df["m_z"]
yinit = df["mobility"]
y1 = df["atDrift1Start"]
y2 = df["atDrift1Stop"]

fit, fitdat, x, y = get_fit(xinit, yinit)
fit1, fitdat1, x1, y1 = get_fit(xinit, y1)
fit2, fitdat2, x2, y2 = get_fit(xinit, y2)

if True:
    calstring = gen_cal(400, fit, fit1, fit2) + gen_cal(900, fit, fit1, fit2) + gen_cal(1400, fit, fit1, fit2)
    write_sonar(calstring)

plt.scatter(x, y, color="k")
plt.plot(fitdat[:, 0], fitdat[:, 1], color="k")
plt.scatter(x1, y1, color="g")
plt.plot(fitdat1[:, 0], fitdat1[:, 1], color="g")
plt.scatter(x2, y2, color="b")
plt.plot(fitdat2[:, 0], fitdat2[:, 1], color="b")
plt.show()
