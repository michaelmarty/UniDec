__author__ = 'michael.marty'

import numpy as np
import matplotlib.pyplot as plt
import json
import os

dir="C:\\cprog\\Me\\Orbitrap\\260115"
os.chdir(dir)

fnames=["CLIP_CA_NEG.json","HIS_CA_NEG.json"]
fnames=["CLIP_NOCA_NEG.json","HIS_NOCA_NEG.json"]
fnames=["CLIP_CA_NEG.json","CLIP_NOCA_NEG.json"]
#fnames=["CLIP_CA_NEG_3.json","HIS_CA_NEG_3.json","CLIP_NOCA_NEG_3.jsdon","HIS_NOCA_NEG_3.json"]
#fnames=["CLIP_CA_NEG_2.json","HIS_CA_NEG_2.json","CLIP_NOCA_NEG_2.json","HIS_NOCA_NEG_2.json"]
#fnames=["CLIP_CA_NEG.json","HIS_CA_NEG.json","CLIP_NOCA_NEG.json","HIS_NOCA_NEG.json"]
#fnames=["collection3.json","collection4.json"]
labels=["POPS","POPC"]
if len(fnames)>2:
    label2=["Ca Clip","Ca His","Clip","His"]
else:
    label2=["No Ca","Ca"]
color=["r","b"]
markers=["o","s",">","v"]
style=["-","--","-.",":"]

plt.figure()
for i in range(0,len(fnames)):
    fname=fnames[i]
    with open(fname,"r") as outfile:
        indict=json.load(outfile)

    extract=np.array(indict["extract"]).transpose()
    print(extract.shape)

    yin=indict["y"]

    yvals=[y[1] for y in yin]

    for j in range(0,len(extract)):
        plt.plot(yvals,extract[j],color=color[j],label=label2[i]+" "+labels[j],marker=markers[i],linestyle=style[i])

    print(yvals)
plt.legend(loc=7)
plt.xlabel("Time (hrs)")
plt.ylabel("Relative Intensity")
#plt.title("No Ca Added")
plt.show()