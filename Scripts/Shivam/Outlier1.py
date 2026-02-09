import PublicScripts.Lipids.LipidOutlierTestScript as ot
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_Dot_prdt_900_RT_Tolerance_0_5_r1_PeakMaster_3_2026_01_29_15_21_48.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_Dot_prdt_900_RT_Tolerance_0_5_r1_PeakMaster_4_2026_01_29_15_16_58.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_Dot_prdt_900_RT_Tolerance_0_5_r1_SpectrumAndReference_AlignmentResult_2026_01_14_10_23_24_2026_01_29_15_21_48.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_Dot_prdt_900_RT_Tolerance_0_5_r_1_SpectrumAndReference_AlignmentResult_2026_01_29_10_01_33_2026_01_29_15_16_58.xml"
os.chdir(topdir)


outdf = ot.outlier_pipeline(posfile, negfile, posxml, negxml, tol=0.3)