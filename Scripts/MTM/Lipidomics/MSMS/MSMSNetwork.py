import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from PublicScripts.Lipids.LipidFunctions import network_plot, network_analysis

mpl.use("WxAgg")

os.chdir(r"Z:\Group Share\Annika\IQ-X\UntDDA Tests HEK\New Libs Tests\Peak Masters")
posfile = "PosPeakMaster_HEKpos_OTIT_WithCorr.csv"
negfile = "NegPeakMaster_HEKneg_OTIT_WithCorr.csv"
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
# posdf.to_csv("PosHEK_OTIT_HeadgroupTailCheck_Network.csv", index=False)
# negdf.to_csv("NegHEK_OTIT_HeadgroupTailCheck_Network.csv", index=False)

# Get global min and max degree for color normalization
colortarget = "z_Anomaly_Score"
colortarget = "Class"

for G in graphs:
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    network_plot(G, ax, colortarget=colortarget)
plt.show()

