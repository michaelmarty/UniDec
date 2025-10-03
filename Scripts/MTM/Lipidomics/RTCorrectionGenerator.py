import pandas as pd
import numpy as np
import os

data = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\NegIDs.csv"
os.chdir("Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX")

df = pd.read_csv(data)

# Need df with "Name, RT (min), RT tolerance (min), m/z (Da), m/z tolerance (Da), Minimum peak height, Include(True or False)"
name = df['Metabolite name']
rt = df['Average Rt(min)']
rtref = df['Reference RT']
mz = df['Average Mz']

# Create new DataFrame
new_df = pd.DataFrame({
    'Name': name,
    'RT (min)': rtref,
    'RT tolerance (min)': 4,
    'm/z (Da)': mz,
    'm/z tolerance (Da)': 0.01,
    'Minimum peak height': 1000,
    'Include(True or False)': True
})

# Save to Text file tab delimited
new_df.to_csv("RTCorrectionLibrary.txt", sep='\t', index=False)