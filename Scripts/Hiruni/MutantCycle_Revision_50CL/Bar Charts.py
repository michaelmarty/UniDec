import pandas as pd
import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams

matplotlib.use('WxAgg')
import matplotlib.pyplot as plt

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.75
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 14
rcParams['font.sans-serif'] = "Arial"



import pandas as pd
import matplotlib.pyplot as plt

#[GA code is here but deleted for simplicity]

x = ["200", "100","50"]
y = ['0.71', '1.40', '1.11']

plt.bar(x,y)

plt.xlabel("CL Concnetration")
plt.ylabel("ΔΔG (kJ/mol)")
plt.title("WT Vs W14A")
plt.show()