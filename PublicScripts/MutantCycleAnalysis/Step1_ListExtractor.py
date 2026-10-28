# Script to extract the peak areas from each deconvolved native mass spectrum
# Note, you must first run UniDec on all the raw data to generate the deconvolved spectra
# This can be done either individually or in batch

import numpy as np
import os
import glob
from unidec.modules.Extract2D import *

R = 8.31446261815324 / 1000.
units = "kJ/mol"

#SpecifyingDirectory
globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final2\\"

os.chdir(globaldir)

#ListExtractorProcessingParameters
edges = [("WT", "W14A")]
temps = [15, 20, 25, 30, 35]
width = 95
method = 6
params1 = [98590, 1400, 350, 0, 7, 0, 1, width, method]
params2 = [98260, 1400, 680, 0, 7, 0, 1, width, method]
params3 = [98260, 1400, 320, 0, 7, 0, 1, width, method]
params4 = [98415, 1400, 500, 0, 7, 0, 1, width, method]
params5 = [98706, 1400, 208, 0, 7, 0, 1, width, method]
params = [params4]

#SpecifyingPath
app = wx.App(False)
for ed in range(0, 1):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All"
    pathtop = os.path.join(globaldir, dirname)
    print(os.path.isdir(pathtop), pathtop)
    out = True

    files = []
    paths = []
    datalist = []
    # Loop through temperatures
    for t in temps:
        newpath = os.path.join(pathtop, e1 + "_" + e2 + "_" + str(t) + "C")
        os.chdir(newpath)
        temp = t
        for path in glob.glob(f'{newpath}/*/**/', recursive=True):
            mass_file = [f for f in os.listdir(path) if "_mass.txt" in f]
            if len(mass_file) > 0:
                mass_file = os.path.join(path, mass_file[0])
                files.append(mass_file)
                paths.append(path)
                data = np.loadtxt(mass_file)
                datalist.append(data)
        print(files)
        if len(files) < 0:
            print("Failed at t=", t)
            break
    print(len(datalist))
    # Load all data into the Extraction GUI, it will automatically export the required text files for Step 2 while
    # Allowing users to check with the window that everything looks good. Just click next through all.
    frame = Extract2DPlot(None, datalist, params=params[ed], directory=paths)
    frame.on_all_next()
    app.MainLoop()
