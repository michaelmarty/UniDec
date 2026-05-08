import os
import pandas as pd
import numpy as np
# d

def tl_to_msp(tl_file, msp_file=None, namecol="Molecule Name", listcol="Molecule List Name",
              formulacol="Molecule Formula", precmzcol="Precursor Mz", rtcol="Molecule Retention Time",
              chargecol="Precursor Charge", adductcol="Precursor Adduct", prodmzcol="Product Mz" ):

    if msp_file is None:
        msp_file = os.path.splitext(tl_file)[0] + ".msp"

    tl_df = pd.read_csv(tl_file)
    if rtcol not in tl_df.columns:
        rtcol = "Explicit Retention Time"
        if rtcol not in tl_df.columns:
            print(f"Warning: Retention time column not found. Setting retention time to 0.")
            tl_df[rtcol] = 0

    # Get unique entries based on the important columns
    udf = tl_df[[namecol, adductcol]].drop_duplicates()

    with open(msp_file, "w") as f:
        for i, row in udf.iterrows():
            # Get the corresponding row in the original dataframe to access all columns
            subdf = tl_df[(tl_df[namecol] == row[namecol]) & (tl_df[adductcol] == row[adductcol])]

            peaks = subdf[prodmzcol].dropna().astype(float).tolist()
            peaks = np.unique(peaks)

            name = row[namecol]
            lclass = subdf[listcol].iloc[0]
            formula = subdf[formulacol].iloc[0]
            precursor_mz = subdf[precmzcol].iloc[0]
            # Convert retention time to seconds if present
            rt_val = np.mean(subdf[rtcol])
            charge = subdf[chargecol].iloc[0]
            if "-" in str(charge):
                polarity = "Negative"
            # elif "+" in str(charge):
            #     polarity = "Positive"
            else:
                polarity = "Positive"
            adduct = row[adductcol]

            f.write(f"NAME: {name}\n")
            f.write(f"PRECURSORMZ: {precursor_mz}\n")
            f.write(f"PRECURSORTYPE: {adduct}\n")
            f.write(f"FORMULA: {formula}\n")
            f.write(f"RETENTIONTIME: {rt_val}\n")
            f.write(f"IONMODE: {polarity}\n")
            f.write(f"COMPOUNDCLASS: {lclass}\n")
            f.write(f"Comment: from {tl_file}\n")
            f.write(f"Num Peaks: {len(peaks)}\n")
            for peak in peaks:
                intensity = 100  # Placeholder intensity, can be modified to use actual intensity if available
                f.write(f"{peak}\t{intensity}\n")
            f.write("\n")  # Separate entries with a blank line
    # Print the file here for confirmation
    print(f"MSP file created: {msp_file}")
    # Print the file contents for confirmation
    with open(msp_file, "r") as f:
        contents = f.read()
        print(contents)


if __name__ == "__main__":
    tl_file = r"Z:\Group Share\Annika\Stellar\FAM\CholTransList.csv"
    tl_file = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData\CleanedAdducts_AllDatasets_TransList_MM_Clean10.csv"
    tl_to_msp(tl_file)
