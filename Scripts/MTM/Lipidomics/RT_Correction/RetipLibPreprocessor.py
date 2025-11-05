import os
from RetipFormatter import *
try:
    import retip
except:
    # add the parent directory to the path to load the Retip library locally in case it isn't installed
    import os, sys
    directory = "C:\\Python\\pyRetip"
    sys.path.insert(1, directory)

    import retip

poslib = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"
posdfname = "Pos_C1_Ref_FromMSP.csv"
neglib = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"
negdfname = "Neg_C1_Ref_FromMSP.csv"

if False:
    os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")
    posdf = msp_to_df(poslib)
    # Write to CSV
    posdf.to_csv(posdfname, index=False)

    # Read neg file
    negdf = msp_to_df(neglib)
    # Write to CSV
    negdf.to_csv(negdfname, index=False)

if True:
    os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
    allowed_columns = ["nAcid", "nBase", "SLogP", "SlogP_VSA5", "SssCH2", "ZMIC1", "CIC0", "MDEC-22", "AATS2pe", "AMID", "CIC1",
     "SMR_VSA5"]

    dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset(testing=posdfname)
    # Limit dataset to first 100 entries for testing
    dataset.calculate_descriptors()
    dataset.describe()
    dataset.preprocess_features('lipidomics')
    dataset.save_retip_dataset(posdfname.replace(".csv", "_preprocessed"))
    dataset.describe()

    dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset(testing=negdfname)
    dataset.calculate_descriptors()
    dataset.describe()
    dataset.preprocess_features('lipidomics')
    dataset.save_retip_dataset(negdfname.replace(".csv", "_preprocessed"))
    dataset.describe()


# Test the df to msp function
# posdf = pd.read_csv("Pos_C1_Ref_FromMSP.csv")
# df_to_msp(posdf, "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp", "Test_Pos_C1.msp")

exit()