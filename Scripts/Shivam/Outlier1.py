import PublicScripts.Lipids.LipidOutlierTestScript as ot
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_updated_PeakMaster_7_2026_01_23_14_30_21.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_updated_PeakMaster_1_2026_01_23_14_07_18.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Positive_Ion_mode\Positive_updated_SpectrumAndReference_AlignmentResult_2026_01_14_17_58_45_2026_01_23_14_30_22.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Negative_Ion_Mode\Negative_updated_SpectrumAndReference_AlignmentResult_2026_01_21_15_54_37_2026_01_23_14_07_18.xml"

os.chdir(topdir)


outdf = ot.outlier_pipeline(posfile, negfile, posxml, negxml, tol=0.3)