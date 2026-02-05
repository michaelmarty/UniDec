import os
import pandas as pd
import numpy as np

# filepath = r"C:\Python\pyRetip\notebooks\lipidomics\PosFullMSMS.csv"
#
# df = pd.read_csv(filepath)
# # Sort df by Name and Adduct to bring duplicates together
# df = df.sort_values(by=['Name', 'Adduct']).reset_index(drop=True)
#
# # Create column to mark duplicates
# df['IsDuplicate'] = df.duplicated(subset=['Name', 'Adduct'], keep=False)
# # Print number of duplicates found
# num_duplicates = df['IsDuplicate'].sum()
# print(f"Number of duplicate entries found: {num_duplicates}")
#
# # Print unique classes that are duplicates
# duplicate_classes = df[df['IsDuplicate']]['CompoundClass'].unique()
# print("Duplicate Classes Found:")
# for cls in duplicate_classes:
#     print(cls)
# exit()
# Write to new excel file
# output_filepath = r"C:\Python\pyRetip\notebooks\lipidomics\NegFullMSMS_DuplicatesMarked.csv"
# df.to_csv(output_filepath, index=False)

filepath = r"C:\Python\pyRetip\notebooks\lipidomics\MSDIAL-TandemMassSpectralAtlas-VS69-Neg-RTP2_Fixed.msp"
filepath = r"C:\Python\pyRetip\notebooks\lipidomics\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-RTP2.msp"
outpath = filepath.replace(".msp", "_Unique.msp")
# Read MSP file
with open(filepath, 'r') as file:
    lines = file.readlines()
    # If line starts with Name add random number to end to make unique
    index = 0
    with open(outpath, 'w') as outfile:
        for line in lines:
            if line.startswith("NAME:"):
                line = line.strip()
                line += f"${index}\n"
                index += 1
            outfile.write(line)