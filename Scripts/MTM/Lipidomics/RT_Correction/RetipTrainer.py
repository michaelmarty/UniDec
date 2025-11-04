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
      'C1_Ref_preprocessed_training.csv',
      'C1_Ref_preprocessed_testing.csv')

drop_cols = ["Reference RT"]
for d in dataset.datasets.keys():
    for col in drop_cols:
        if col in dataset.datasets[d].columns:
            dataset.datasets[d] = dataset.datasets[d].drop(columns=[col])

dataset.describe()
# Print dataset columns
print(dataset.datasets["training"].columns)
dataset.descriptor_names.append('S/N average')
dataset.descriptor_names.append("Tail Length")
dataset.descriptor_names.append("Tail Unsaturation")


# trainer = retip.XGBoostTrainer(dataset)
trainer = retip.AutoGluonTrainer(dataset, training_duration=3, preset="high_quality")
trainer.weight_column = "S/N average"
trainer.train()
print(trainer.model_columns)

trainer.save_model('C1_Ref_autogluon-model_HighW10')

print(trainer.score(plot=False))