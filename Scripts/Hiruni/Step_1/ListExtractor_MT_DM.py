import numpy as np
import os

topdir = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final\\R75A_R224-75A Titrations All\\R224A_R224-75A_40C\\"

files = ["rep1\\202200510_HSJ_AqpZR75A_R224-75A_100CL_TS_40_1_unidecfiles", "rep2\\202200510_HSJ_AqpZR75A_R224-75A_100CL_TS_40_2_unidecfiles", "rep3\\202200510_HSJ_AqpZR75A_R224-75A_100CL_TS_40_3_unidecfiles"]

datalist = []
paths = []
for f in files:
    path = os.path.join(topdir, f)
    paths.append(path)
    os.chdir(path)
    print(path, f)
    mass_file= [f for f in os.listdir(path) if "_mass.txt" in f][0]
    print (mass_file)
    data = np.loadtxt(mass_file)
    datalist.append(data)

app = wx.App(False)
params = [98260, 1400, 300, 0, 7, 0, 1, 100, 2]
frame = Extract2DPlot(None, datalist, params=params, directory=paths)

# frame.run_multip("1,2,3,4,5,6")
app.MainLoop()