import numpy as np
import os
import glob
from unidec.modules.Extract2D import *

R = 8.31446261815324 / 1000.
units = "kJ/mol"

globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\WT_K4A\\"

os.chdir(globaldir)
edges = [("WT", "K4A")]
temps = [15, 20, 25, 30, 35, 40]
width = 70
method = 6
params1 = [98590, 1400, 350, 0, 7, 0, 1, width, method]
params2 = [98260, 1400, 680, 0, 7, 0, 1, width, method]
params3 = [98260, 1400, 300, 0, 7, 0, 1, width, method]
params4 = [98400, 1400, 500, 0, 7, 0, 1, width, method]
params5 = [98700, 1400, 205, 0, 7, 0, 1, width, method]
params = [params5]
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
        if len(files) < 1:
            print("Failed at t=", t)
            break
    print(len(datalist))

    frame = Extract2DPlot(None, datalist, params=params[ed], directory=paths)
    frame.on_all_next()
    app.MainLoop()
