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


def split_csvs(inputcsv, outputpath):
    df = pd.read_csv(inputcsv)
    classes = df['CompoundClass'].unique()
    outfiles = []
    for c in classes:
        class_df = df[df['CompoundClass'] == c]
        class_df.to_csv(os.path.join(outputpath, f"{c.replace('/', '_')}_data.csv"), index=False)
        outfiles.append(os.path.join(outputpath, f"{c.replace('/', '_')}_data.csv"))
    return outfiles


def run_preprocessor(input_csv, output_csv, allowed_columns=None):
    if os.path.isfile(output_csv + "_training.csv"):
        print(f"Output file {output_csv}_training.csv already exists. Skipping preprocessing.")
        return
    dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset(input_csv)
    dataset.calculate_descriptors(n_proc=12)
    dataset.save_retip_dataset(output_csv)
    dataset.describe()

def run_preprocessor_batch(input_csv, output_csv, allowed_columns=None, batchsize=1000):
    if os.path.isfile(output_csv + "_training.csv"):
        print(f"Output file {output_csv}_training.csv already exists. Skipping preprocessing.")
        return
    dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset(input_csv)

    topdf = dataset.datasets["training"].copy()
    print(topdf)

    for i in range(0, len(topdf), batchsize):
        batch_output_csv = output_csv.replace("_preprocessed_full", f"_batch{i // batchsize + 1}")
        # Check if file exists
        if os.path.isfile(batch_output_csv + "_training.csv"):
            print(f"Batch output file {batch_output_csv}_training.csv already exists. Skipping this batch.")
            continue

        batch_df = topdf.iloc[i:i+batchsize]
        temp = batch_df.copy()
        temp.to_csv("temp_batch.csv", index=False)
        batch_dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset("temp_batch.csv")
        # batch_dataset.datasets["training"] = batch_df
        batch_dataset.calculate_descriptors(n_proc=20)
        batch_dataset.save_retip_dataset(batch_output_csv)
        batch_dataset.describe()

    # Stitch them back together
    for i in range(0, len(topdf), batchsize):
        batch_output_csv = output_csv.replace("_preprocessed_full", f"_batch{i // batchsize + 1}") + "_training.csv"
        batch_dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns).load_retip_dataset(batch_output_csv)
        if i == 0:
            final_df = batch_dataset.datasets["training"]
        else:
            final_df = pd.concat([final_df, batch_dataset.datasets["training"]], ignore_index=True)

    final_dataset = retip.Dataset(target_column='RT', allowed_descriptors=allowed_columns)
    final_dataset.datasets["training"] = final_df
    final_dataset.save_retip_dataset(output_csv)
    final_dataset.describe()

posdfname = "Pos_C1_Ref_FromMSP.csv"
negdfname = "Neg_C1_Ref_FromMSP.csv"

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics\SplitCSVs")

file = "Pos/Vitamin_D_data.csv"
outputfile = file.replace(".csv", "_preprocessed_full")
print(f"Processing {file} -> {outputfile}")
run_preprocessor_batch(file, outputfile)




