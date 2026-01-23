import PublicScripts.Lipids.LipidOutlierTestScript as ot
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks"
posfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Whole_class_of_lipids_positive_ion_mode.csv"
negfile = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Neg_ion_using_fixed_libraryPeakMaster_1_2026_01_21_16_05_03.csv"
posxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Pos_Ion_mode_SpectrumAndReference_AlignmentResult_2026_01_14_17_58_45_2026_01_20_12_57_28.xml"
negxml = r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks\Neg_ion_using_fixed_library_SpectrumAndReference_AlignmentResult_2026_01_21_15_54_37_2026_01_21_16_05_14.xml"

os.chdir(topdir)

outdf = ot.outlier_pipeline(posfile, negfile, posxml, negxml, tol=0.3)