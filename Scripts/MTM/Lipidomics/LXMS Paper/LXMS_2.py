import os
from scipy import stats
import pingouin as pg

if True:
    os.chdir("C:\Data\Lipidomics\\AchR")

    # Set Paths
    path = "LXMS_AchR_Output.xlsx"

    outfile = "LXMS_AchR_Classes.xlsx"

    drop_classes = ["LPC", "LPE", "LPE-O", "LPE-O/P", "TG", "PG", "Cer", "CL"]
    sample_types = ["FT_Ctrl", "Elu_Ctrl", "FT_AchR", "Elu_AchR", "E3(+)_Ctrl", "E3(+)_AchR"]
    outliers = ["FT_Ctrl_1"]
    outliers = []

    comparisons = [["FT_Ctrl", "Elu_Ctrl"],
                   ["E3(+)_Ctrl", "E3(+)_AchR"],
                   ["FT_AchR", "Elu_AchR"],
                   ]

else:
    os.chdir("C:\Data\Lipidomics\\AmtBAqpZ")
    path = "LXMS_Ecoli_Output.xlsx"

    outfile = "LXMS_Ecoli_Classes.xlsx"

    drop_classes = ["LPE", "LPE-O", "LPE-O/P", ]

    sample_types = ["Elu_Ctrl", "FT_Ctrl", "Elu_AmtB", "FT_AmtB", "E3_Ctrl", "E3_AmtB",
                    "Elu_AqpZ", "FT_AqpZ", "E3_Ctrl", "E3_AqpZ"]
    outliers = []

    comparisons = [["Elu_Ctrl", "FT_Ctrl"],
                   ["E3_Ctrl", "E3_AmtB"],
                   ["Elu_AmtB", "FT_AmtB"],
                   ["E3_Ctrl", "E3_AqpZ"],
                   ["Elu_AqpZ", "FT_AqpZ"],
                   ]
# sample_types = ["4_D1_Ecoli", "4_E3_Ecoli", "25_D1_Ecoli", "25_E3_Ecoli", "37_D1_Ecoli", "37_E3_Ecoli"]


sample_numers = ["1", "2", "3"]

# Load the data
df = pd.read_excel(path)

# Strip off the class from the name
classes = df["Class Name"]

# Drop all -O/P from class names
classes2 = []
for c in classes:
    if "IS" not in c:
        c = re.sub(r'-O/P', '', c)
    classes2.append(c)
classes = np.array(classes2)

uclasses = np.unique(classes)

all_classes = np.unique(classes)

print(uclasses)

# Drop classes from the uclasses array

# drop_classes = []
for c in uclasses:
    if "IS" in c:
        drop_classes.append(c)

for drop_class in drop_classes:
    uclasses = uclasses[uclasses != drop_class]

# Move SM in the uclasses array to after HexCer
if "SM" in uclasses and "HexCer" in uclasses:
    insert_index = np.where(uclasses == "HexCer")[0][0] + 1
    uclasses = np.delete(uclasses, np.where(uclasses == "SM"))
    uclasses = np.insert(uclasses, insert_index, "SM")
# Move TG to the front
if "TG" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "TG"))
    uclasses = np.insert(uclasses, 0, "TG")

# Move CL to the back
if "CL" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "CL"))
    uclasses = np.insert(uclasses, len(uclasses), "CL")

# include_classes = ["PE-O"]
# uclasses = include_classes
colors = get_colors(uclasses)

# All samples
samples = []
class_sums = []
class_tls = []
class_tus = []

# Sample means and stds
class_means = []
class_stds = []
tls_means = []
tls_stds = []
tus_means = []
tus_stds = []

if "Pos" in df and "Neg" in df:
    check_mode = True
    posmode = df["Pos"]
else:
    check_mode = False
    posmode = np.ones_like(df["Class Name"])

for sample_type in sample_types:

    sample_classes = []
    sample_tls = []
    sample_tus = []

    for sample_num in sample_numers:

        sample_name = sample_type + "_" + sample_num

        if sample_name in outliers:
            continue

        print(sample_name)
        samples.append(sample_name)

        # Extract data frame where column names are in exp_name
        exp_df = df[sample_name]

        # exp_df = df[df.columns.intersection([sample_name])]

        # Set a boolean of all classes that are in uclasses
        class_bool_smpl = np.array([c in uclasses for c in classes])
        # Set a boolean of all classes that are internal standards
        class_is_sample = np.array(["IS" in c for c in classes])

        exp_mtic = float(np.sum(exp_df[class_bool_smpl]))

        exp_stic = float(np.sum(exp_df[class_is_sample]))

        exp_normed_sum = []
        tl_class = []
        tu_class = []

        for uclass in uclasses:
            # Get main class bool
            class_bool = classes == uclass
            # Get internal standard class bool
            isclass = get_internal_standards(uclass, all_classes)
            class_is = classes == isclass
            # Calculate total internal standard intensity
            exp_df_is = exp_df[class_is]

            exp_abs = exp_df[class_bool].to_numpy()  # / exp_mtic

            if len(exp_df_is) > 2:
                print("Error in internal standard, more than two found", isclass, uclass, all_classes)
                raise Exception
            elif len(exp_df_is) == 2:

                exp_pol = posmode[class_bool].to_numpy()
                is_pol = posmode[class_is].to_numpy()

                exp_is = np.ones_like(exp_abs)
                exp_is[exp_pol == 0] = exp_df_is[is_pol == 0]
                exp_is[exp_pol == 1] = exp_df_is[is_pol == 1]

            elif len(exp_df_is) == 1:
                if check_mode:
                    exp_pol = posmode[class_bool].to_numpy()
                    is_pol = posmode[class_is].to_numpy()[0]
                    if np.any(exp_pol != is_pol):
                        pol_match_bool = exp_pol == is_pol
                        # Remove non matching polarities
                        exp_abs = exp_abs[pol_match_bool]
                        class_bool[class_bool] *= pol_match_bool

                        print("Error in internal standard, polarity does not match", isclass, uclass)
                        # raise Exception

                exp_is = float(np.sum(exp_df_is))  # / exp_stic
            # Normalize by internal standard
            if np.sum(exp_is) == 0:
                print("No internal standard for class: ", uclass)
                raise Exception

            exp_normed = exp_abs / exp_is

            exp_normed_sum.append(np.sum(exp_normed))

            cdf = df[class_bool]

            tl = cdf["Tail Lengths"]
            tu = cdf["Tail Unsaturation"]
            tld = cdf["Tail Length Diff"]
            tud = cdf["Tail Unsat Diff"]

            bool_tl = tl >= 0
            bool_tu = tu >= 0

            exp_avg_tl = np.average(tl[bool_tl], weights=np.sqrt(exp_normed[bool_tl]))
            exp_avg_tu = np.average(tu[bool_tu], weights=np.sqrt(exp_normed[bool_tu]))

            tl_class.append(exp_avg_tl)
            tu_class.append(exp_avg_tu)

        total_exp = np.sum(exp_normed_sum)
        exp_normed_sum = np.array(exp_normed_sum)
        exp_normed_sum_percent = exp_normed_sum / total_exp

        class_sums.append(exp_normed_sum_percent)
        class_tls.append(tl_class)
        class_tus.append(tu_class)

        sample_classes.append(exp_normed_sum_percent)
        sample_tls.append(tl_class)
        sample_tus.append(tu_class)

    class_means.append(np.mean(sample_classes, axis=0))
    class_stds.append(np.std(sample_classes, axis=0))
    tls_means.append(np.mean(sample_tls, axis=0))
    tls_stds.append(np.std(sample_tls, axis=0))
    tus_means.append(np.mean(sample_tus, axis=0))
    tus_stds.append(np.std(sample_tus, axis=0))

writer = pd.ExcelWriter(outfile, engine='openpyxl')
# Save the results
class_df1 = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df1[sample] = class_sums[i]

class_df1.to_excel(writer, index=False, sheet_name="Classes")

# Save the results
class_df2 = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df2[sample] = class_tls[i]

class_df2.to_excel(writer, index=False, sheet_name="Tail Lengths")

# Save the results
class_df3 = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df3[sample] = class_tus[i]

class_df3.to_excel(writer, index=False, sheet_name="Tail Unsaturation")

# Save the results
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(sample_types):
    class_df[sample + " Mean"] = class_means[i]
    class_df[sample + " Std"] = class_stds[i]

class_df.to_excel(writer, index=False, sheet_name="Class Means")

# Save the results
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(sample_types):
    class_df[sample + " Mean"] = tls_means[i]
    class_df[sample + " Std"] = tls_stds[i]

class_df.to_excel(writer, index=False, sheet_name="Tails Means")

# Save the results
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(sample_types):
    class_df[sample + " Mean"] = tus_means[i]
    class_df[sample + " Std"] = tus_stds[i]

class_df.to_excel(writer, index=False, sheet_name="Unsaturation Means")

output_data = []
for i, col in enumerate(samples):
    splits = col.split("_")
    lx_type = splits[0]
    nd_type = splits[1]
    sample_number = splits[2]
    total_type = lx_type + "_" + nd_type

    print(col, lx_type, nd_type)

    for j, c in enumerate(uclasses):
        intensity = class_sums[i][j]
        tl = class_tls[i][j]
        tu = class_tus[i][j]
        output_data.append([col, total_type, lx_type, nd_type, sample_number, c, intensity, tl, tu])

outdf = pd.DataFrame(output_data,
                     columns=["Sample Name", "Sample Type", "LX_Type", "ND_Type", "Sample Number",
                              "Class Name", "Intensity", "Tail Length", "Tail Unsaturation"])

outdf.to_excel(writer, index=False, sheet_name="Reformat")

rawdfs = [class_df1, class_df2, class_df3]
int_names = ["Intensity", "Tail Length", "Tail Unsaturation"]
cutoffs = [0.05, 0.05, 0.05]
# rawdfs = [class_df1]

for x, rawdf in enumerate(rawdfs):
    results = []
    for i, comparison in enumerate(comparisons):
        print(comparison)

        column_names1 = [col for col in rawdf.columns if comparison[0] in col]
        column_names2 = [col for col in rawdf.columns if comparison[1] in col]

        df1 = rawdf[column_names1]
        df2 = rawdf[column_names2]

        if x == 0:
            df1 = df1 * 100
            df2 = df2 * 100

        m1 = np.mean(df1, axis=1)
        s1 = np.std(df1, axis=1)
        m2 = np.mean(df2, axis=1)
        s2 = np.std(df2, axis=1)

        fc = m2 / m1
        fcsd = np.sqrt((s1 / m1) ** 2 + (s2 / m2) ** 2) * fc

        diff = m2 - m1
        diffsd = np.sqrt(s1 ** 2 + s2 ** 2)

        avg_fc = np.mean(fc)
        std_fc = np.std(fc)

        avg_diff = np.mean(diff)
        std_diff = np.std(diff)

        normality1 = np.array([stats.shapiro(df1.iloc[i])[1] for i in range(len(df1))])
        normality2 = np.array([stats.shapiro(df2.iloc[i])[1] for i in range(len(df2))])

        b1 = normality1 < 0.05
        b2 = normality2 < 0.05

        print("Normality 1: ", normality1[b1], "Normality 2: ", normality2[b2])

        pvals = stats.ttest_ind(df1, df2, axis=1)[1]

        for j, c in enumerate(uclasses):
            results.append([comparison[1], comparison[0], c, m1[j], s1[j], m2[j], s2[j],
                            fc[j], fcsd[j], diff[j], diffsd[j], pvals[j]])

    results = np.array(results)
    pvals_all = results[:, -1].astype(float)

    if x == 0:
        fc_all = results[:, -5].astype(float)
        b1 = np.abs(fc_all - 1) < cutoffs[x]
        pvals_all[b1] = 2
    else:
        diffs_all = results[:, -3].astype(float)
        b1 = np.abs(diffs_all) < cutoffs[x]
        pvals_all[b1] = 2
    # b2 is not b1
    b2 = pvals_all != 2
    print(np.sum(b1), np.sum(b2))
    adjusted_pvals = pvals_all
    adjusted_pvals[b2] = pg.multicomp(pvals_all[b2], method="fdr_bh")[1]

    full_results = pd.DataFrame(results, columns=["A", "B", "Class Name", "A_Mean", "A_Std", "B_Mean", "B_Std",
                                                  "FC", "FC_Std", "Diff", "Diff_Std", "pval"])
    full_results["Adjusted pval"] = adjusted_pvals

    full_results.to_excel(writer, index=False, sheet_name=int_names[x] + " Stats")

writer.close()
