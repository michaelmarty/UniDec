import PublicScripts.Lipids.LipidFunctions as lf
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Whole_class_of_lipids_positive_ion_mode.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Whole_class_lipids_PeakMaster_neg_ion_mode.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Pos_Ion_mode_SpectrumAndReference_AlignmentResult_2026_01_14_17_58_45_2026_01_20_12_57_28.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Neg_Ion_mode_SpectrumAndReference_AlignmentResult_2026_01_16_11_08_43_2026_01_20_13_01_23.xml"


mpl.use("WxAgg")
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

posdf, posgraphs = lf.class_network_analysis(posdf)
negdf, neggraphs = lf.class_network_analysis(negdf)
graphs = posgraphs + neggraphs

# posdf, posgraphs = lf.network_analysis(posdf)
# negdf, neggraphs = lf.network_analysis(negdf)
# graphs = [neggraphs] + [posgraphs]

# posdf, posgraphs = lf.tail_network_analysis(posdf, min_count=5)
# negdf, neggraphs = lf.tail_network_analysis(negdf, min_count=5)
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
    lf.network_plot(G, ax, colortarget=colortarget)
plt.show()