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

trainer = retip.AutoGluonTrainer(dataset)
trainer.load_model('C1_Ref_autogluon-model_HighW10')

annotated = trainer.annotate(dataset.get_testing_data(include_metadata=True), prediction_column='RTP')
# Save
annotated.to_csv("C1_Ref_Testing_HighW10.csv")

print(trainer.score(plot=True))

retip.visualization.plot_feature_importance(trainer)


