import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplcursors

from Scripts.MTM.Lipidomics.LipidFunctions import class_network_analysis, network_plot, network_analysis, tail_network_analysis

mpl.use("WxAgg")

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
posfile = "PosIDsRTP2_HeadgroupTailCheck.csv"
negfile = "NegIDsRTP2_HeadgroupTailCheck.csv"
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

# posdf, posgraphs = class_network_analysis(posdf)
# negdf, neggraphs = class_network_analysis(negdf)
# graphs = posgraphs + neggraphs

posdf, posgraphs = network_analysis(posdf)
negdf, neggraphs = network_analysis(negdf)
graphs = [neggraphs] + [posgraphs]

# posdf, posgraphs = tail_network_analysis(posdf, min_count=5)
# negdf, neggraphs = tail_network_analysis(negdf, min_count=5)
# graphs = posgraphs + neggraphs

# Save updated dataframes
posdf.to_csv("PosIDsRTP2_HeadgroupTailCheck_Network.csv", index=False)
negdf.to_csv("NegIDsRTP2_HeadgroupTailCheck_Network.csv", index=False)

# Get global min and max degree for color normalization
colortarget = "z_Anomaly_Score"
colortarget = "Class"

for G in graphs:
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    network_plot(G, ax, colortarget=colortarget)
plt.show()

