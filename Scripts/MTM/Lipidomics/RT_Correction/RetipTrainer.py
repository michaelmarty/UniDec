import os
import numpy as np

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
      'Filtered_C1_Ref_Testing_HighWLimited30.csv',
      'C1_Ref_preprocessed_testing.csv')

dataset.describe()
# Print dataset columns
print("Columns in Training Data:", dataset.datasets["training"].columns)

# The top 12 features from a more extensive training plus nAcid and nBase, which I think are important
allowed_columns = ["nAcid", "nBase", "SLogP", "SlogP_VSA5", "SssCH2", "ZMIC1", "CIC0", "MDEC-22", "AATS2pe", "AMID", "CIC1", "SMR_VSA5"]

new_descriptor_names = []
for col in allowed_columns:
    if col in dataset.descriptor_names:
        new_descriptor_names.append(col)

# new_descriptor_names.append('S/N average')
# new_descriptor_names.append("Tail Length")
# new_descriptor_names.append("Tail Unsaturation")
dataset.descriptor_names = new_descriptor_names

# Use log of weight column
# dataset.datasets["training"]["S/N average"] = dataset.datasets["training"]["S/N average"].apply(lambda x: np.log2(x))

trainer = retip.AutoGluonTrainer(dataset, training_duration=30, preset="medium_quality")
# trainer.weight_column = "S/N average"
trainer.train()
print(trainer.model_columns)

trainer.save_model('C1_Ref_autogluon-model_FilteredMedLimited30')
#
# print(trainer.score(plot=False))