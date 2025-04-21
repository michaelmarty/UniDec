import numpy as np
import os
import glob
from unidec.modules.Extract2D import *

R = 8.31446261815324 / 1000.
units = "kJ/mol"

globaldir = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_MAE_MutantCyclePaper\\HSJ_MAE AqpZ CL\\MCP_Revisions_05202024\\Incubation Times\\"

os.chdir(globaldir)
edges = [("WT", "W206A")]
temps = [20, 25, 30]
width = 80
method = 6
params1 = [98425, 1410, 497, 0, 7, 0, 1, width, method] #WT_W14A
params2 = [98560, 1405, 340, 0, 7, 0, 1, width, method] #WT_R224A/R75A
params8 = [98405, 1405, 160, 0, 7, 0, 1, width, method] #R224A_W14A
params = [params1]
app = wx.App(False)
for ed in range(0, 1):
    e1 = edges[ed][0]
    e2 = edges[ed][1]
    #if e2 is "DM":
        #e2 = "R224-75A"
    dirname = e1 + "_" + e2 + " Titrations All_ON"
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