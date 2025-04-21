import pandas as pd
import numpy as np

samples = [1,2,3,4,5,6,7,8,9]
reps = [1,2,3]

rep1 = pd.read_excel("Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\Rep1\\Normalized_Abundances_Rep1.xlsx")
rep2 = pd.read_excel("Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\Rep2\\Normalized_Abundances_Rep2.xlsx")
rep3 = pd.read_excel("Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\Rep3\\Normalized_Abundances_Rep3.xlsx")

names1 = rep1["Lipid Name"].to_numpy()
names2 = rep2["Lipid Name"].to_numpy()
names3 = rep3["Lipid Name"].to_numpy()



names = np.intersect1d(np.intersect1d(names1, names2), names3)

output_df = pd.DataFrame()

output_df["Lipid Name"] = names

#initialize the output dataframe columns
for i in range(len(samples)):
    for j in range(len(reps)):
        output_df["S" + str(samples[i]) + "_R" + str(reps[j])] = np.zeros(len(names))
    output_df["S" + str(samples[i]) + "_average"] = np.zeros(len(names))
    output_df["S" + str(samples[i]) + "_StdDev"] = np.zeros(len(names))


for ni in range(len(names)):
    rep1_index = np.where(rep1["Lipid Name"] == names[ni])[0][0]
    rep2_index = np.where(rep2["Lipid Name"] == names[ni])[0][0]
    rep3_index = np.where(rep3["Lipid Name"] == names[ni])[0][0]
    for i in range(len(samples)):
        areas = []
        for j in range(len(reps)):
            if reps[j] == 1:
                output_df.loc[ni, "S" + str(samples[i]) + "_R" + str(reps[j])] = rep1[str(samples[i])][rep1_index]
                areas.append(rep1[str(samples[i])][rep1_index])
            elif reps[j] == 2:
                output_df.loc[ni, "S" + str(samples[i]) + "_R" + str(reps[j])] = rep2[str(samples[i])][rep2_index]
                areas.append(rep2[str(samples[i])][rep2_index])
            elif reps[j] == 3:
                output_df.loc[ni, "S" + str(samples[i]) + "_R" + str(reps[j])] = rep3[str(samples[i])][rep3_index]
                areas.append(rep3[str(samples[i])][rep3_index])
        output_df.loc[ni, "S" + str(samples[i]) + "_average"] = np.mean(areas)
        output_df.loc[ni, "S" + str(samples[i]) + "_StdDev"] = np.std(areas)

output_df.to_excel('D:\\Johnny\\Lipidomics\\Collaborator_Samples\\Averaged_Normalized_Abundances.xlsx', index=False)

