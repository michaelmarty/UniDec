import PublicScripts.Lipids.LipidFunctions as lf
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_Dot_prdt_900_RT_Tolerance_0_5_r1_PeakMaster_3_2026_01_29_15_21_48.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_Dot_prdt_900_RT_Tolerance_0_5_r1_PeakMaster_4_2026_01_29_15_16_58.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_Dot_prdt_900_RT_Tolerance_0_5_r1_SpectrumAndReference_AlignmentResult_2026_01_14_10_23_24_2026_01_29_15_21_48.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_Dot_prdt_900_RT_Tolerance_0_5_r_1_SpectrumAndReference_AlignmentResult_2026_01_29_10_01_33_2026_01_29_15_16_58.xml"

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