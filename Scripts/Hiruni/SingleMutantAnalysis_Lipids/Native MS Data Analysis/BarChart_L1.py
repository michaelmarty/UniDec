import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from biopandas import pdb
from Bio import SeqIO, SeqUtils
import Bio.PDB.Polypeptide as pp
import os, fnmatch


import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm

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

rcParams['lines.markersize'] = 7
rcParams['font.size'] = 16
rcParams['font.sans-serif'] = "Arial"


def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1))/255
    c2_rgb = np.array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]

path= os.path.abspath("Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis2\\")
datafile = "Z:\\Group Share\\Hiruni Jayasekera\\FAM_SingleMutantAnalysis2\\POPG\\CombinedResults_POCL_Difference.xlsx"

df = pd.read_excel(datafile)
y = df["GA"]

color1 = "#FF0000"
color2 = "#FFFFFF"
plt.gca().set_ylim(-1.0, 1.0)
plt.yticks([-1.0, -0.5, 0, 0.5, 1.0])

#colors = get_color_gradient(color1, color2, len(y))

cmap = plt.get_cmap("RdBu_r")
colors = [cmap((x+1)/2) for x in y]

plt.bar(df["P2"], y, yerr=df["GAE"], color = colors, width = 0.8, linewidth =1, capsize=4)
plt.xticks(rotation=90, ha='right')

plt.savefig(os.path.join(path, "POCL_Difference_0.5.pdf"))

plt.show()

