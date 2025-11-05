import os

try:
    import retip
except:
    # add the parent directory to the path to load the Retip library locally in case it isn't installed
    import os, sys
    directory = "C:\\Python\\pyRetip"
    sys.path.insert(1, directory)

    import retip

import matplotlib.pyplot as plt

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

fname = 'Pos_C1_Ref_FromMSP_preprocessed_testing.csv'

dataset = retip.Dataset(target_column='RT').load_retip_dataset(training='Filtered_C1_Ref_Testing_HighWLimited30.csv',
      testing=fname)

trainer = retip.AutoGluonTrainer(dataset)
trainer.load_model('C1_Ref_autogluon-model_FilteredHighLimited10')

annotated = trainer.annotate(dataset.get_testing_data(include_metadata=True), prediction_column='RTP')
# Save
annotated.to_csv(fname.replace("_preprocessed_testing", "_predicted"))

plt.scatter(annotated['RT'], annotated['RTP'])
plt.xlabel("Experimental RT (min)")
plt.ylabel("Predicted RT (min)")
plt.plot([0, max(annotated['RT'])], [0, max(annotated['RT'])], 'r--')
plt.show()


