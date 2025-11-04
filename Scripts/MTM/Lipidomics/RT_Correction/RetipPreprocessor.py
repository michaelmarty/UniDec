import os

try:
    import retip
except:
    # add the parent directory to the path to load the Retip library locally in case it isn't installed
    import os, sys
    directory = "C:\\Python\\pyRetip"
    sys.path.insert(1, directory)

    import retip

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

dataset = retip.Dataset(target_column='RT').load_retip_dataset(
    training='C1_Ref_Training.csv',
    testing='C1_Ref_Training.csv')

# dataset = retip.Dataset(target_column='RT').load_retip_dataset(
#     training='C1_Ref_Training.csv')
#
# dataset.split_dataset(0.1)

dataset.calculate_descriptors()
dataset.describe()
dataset.preprocess_features('lipidomics')
dataset.save_retip_dataset('C1_Ref_preprocessed')

dataset.describe()

# dataset = retip.Dataset(target_column='RT').load_retip_dataset(
#      'C1_Ref_preprocessed_testing.csv',
#      'C1_Ref_preprocessed_training.csv')