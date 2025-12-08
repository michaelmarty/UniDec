import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def cleanup_df(df):
    # Remove rows with null in MS/MS Spectrum
    df = df[df['MSMS spectrum'].notnull()]
    # Remove name if Unknown
    df = df[~df['Name'].str.contains("Unknown", na=False)]
    df = df[~df['Name'].str.contains("low score", na=False)]
    return df

os.chdir(r"Z:\Group Share\Annika\Stellar\HEK")
file = r"20250929_AS_HEK_DDA_20_Neg.txt"
file = r"20251027_AS_HEK_DDA_Neg_Stellar_RTCor2.txt"
os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")
file = "NegIDsRTP1.csv"

# df = pd.read_table(file)
df = pd.read_csv(file)
# df = cleanup_df(df)

# Plot Reference RT vs Average Rt(min) for each lipid class
classes = df['Ontology'].unique()
plt.figure(figsize=(10, 6))
index = 0
classlist = []
for i, lipid_class in enumerate(classes):
    class_df = df[df['Ontology'] == lipid_class]
    if len(class_df) < 3:
        continue
    index += 1
    classlist.append(lipid_class)
    # plt.scatter(class_df['Reference RT'], class_df['RT (min)'], label=lipid_class)
    # Make a violin plot instead
    diffs = class_df['Average Rt(min)'] - class_df['Reference RT']
    average_rt = np.mean(class_df['Average Rt(min)'])
    plt.subplot(121)
    plt.violinplot(diffs, positions=[index], showmeans=True)
    plt.scatter([index]*len(diffs), diffs, alpha=0.5, s=10)
    plt.subplot(122)
    plt.plot(average_rt, np.mean(diffs), 'o', markersize=len(class_df), label=lipid_class)

plt.subplot(121)
plt.xticks(range(len(classlist)), classlist, rotation=90)

plt.title('RT Comparison by Lipid Class')
plt.show()