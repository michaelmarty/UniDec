import PublicScripts.Lipids.LipidOutlierTestScript as ot
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Whole_class_of_lipids_positive_ion_mode.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Whole_class_lipids_PeakMaster_neg_ion_mode.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Pos_Ion_mode_SpectrumAndReference_AlignmentResult_2026_01_14_17_58_45_2026_01_20_12_57_28.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Neg_Ion_mode_SpectrumAndReference_AlignmentResult_2026_01_16_11_08_43_2026_01_20_13_01_23.xml"

os.chdir(topdir)

outdf = ot.outlier_pipeline(posfile, negfile, posxml, negxml, tol=0.3)