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
    dataset.calculate_descriptors()
    dataset.save_retip_dataset(output_csv)
    dataset.describe()


posdfname = "Pos_C1_Ref_FromMSP.csv"
negdfname = "Neg_C1_Ref_FromMSP.csv"

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics\SplitCSVs")

file = "Pos/TG_data.csv"
outputfile = file.replace(".csv", "_preprocessed_full")
print(f"Processing {file} -> {outputfile}")
run_preprocessor(file, outputfile)




