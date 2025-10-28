import numpy as np
import pandas as pd
import os
from unidec.LipiDec.AddFeaturestoDF import *
from scipy import stats

os.chdir("Z:\\Group Share\\Guozhi Zhang\\2D_LC\\230425_maybe_LX_paper")

# Set Paths
path = "EColi_LXMS_pmartR_output_230426_NEW-mTIC_log2.xlsx"

# Load the data
df = pd.read_excel(path)

df = set_tails(df, "Name")

df = set_class_name(df, "Name")



# Read column names
namedf = pd.read_csv("230426_Ecoli_ND_tag_comp_fdata_1.csv")

comparisons = [["Lipid", "ND_minus"], ["Lipid", "ND_plus"],
               ["ND_minus", "ND_plus"]]

for comp in comparisons:
    label = comp[0] + "-" + comp[1]
    # Get Row where the Label column is Elu
    ctl_name = namedf[namedf["Combined"] == comp[0]]["SampleName"].to_numpy()
    exp_name = namedf[namedf["Combined"] == comp[1]]["SampleName"].to_numpy()

    # Extract data frame where column names are in ctl_name
    ctl_df = df[df.columns.intersection(ctl_name)]
    # Extract data frame where column names are in exp_name
    exp_df = df[df.columns.intersection(exp_name)]

    print(label, len(ctl_df), len(exp_df))
    # Calculate the mean of each df
    ctl_mean = ctl_df.mean(axis=1)
    exp_mean = exp_df.mean(axis=1)

    df["Mean " + label + " 0"] = pd.Series(ctl_mean)
    df["Mean " + label + " 1"] = pd.Series(exp_mean)

    # Calculate the fold change
    fold_change = exp_mean - ctl_mean

    # Add the fold change to the original data frame
    df["Fold Change " + label] = pd.Series(fold_change)

    # Get P value
    p_value = np.zeros(len(df))
    for i, row in df.iterrows():
        p_value[i] = stats.ttest_ind(exp_df.iloc[i, :], ctl_df.iloc[i, :], nan_policy="omit")[1]

    # Add the p value to the original data frame
    df["P Value " + label] = pd.Series(p_value)
    # Print number of nans
    print("Nans: ", np.sum(np.isnan(p_value)))

# labels = ["pre_Ex", "FT", "Elu", "E3(+)"]

# For all columns with Mean in the title, average them together
# Get all columns with Mean in the title
mean_cols = [col for col in df.columns if "Mean" in col]
# Get subdf with only the mean columns
mean_df = df[mean_cols]
# Get the mean of each row
mean_df = mean_df.mean(axis=1)
# Add the mean to the original data frame
df["Mean Abundance Log2"] = pd.Series(mean_df)
# Add the mean after undoing the log2
df["Mean Abundance"] = 2 ** mean_df

df = calc_tail_diff(df, weight_col="Mean Abundance")


# Save the results
df.to_excel("230426_EColi_LXMS_Processed.xlsx", index=False)
