import os

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_Brain_Both_output.xlsx"
path = "NA_Ecoli_output.xlsx"

outfile = "NA_Brain_Classes.xlsx"
outfile = "NA_Ecoli_Classes.xlsx"

sample_types = ["4_D1_Brain", "4_E3_Brain", "25_D1_Brain", "25_E3_Brain", "37_D1_Brain", "37_E3_Brain"]
sample_types = ["4_D1_Ecoli", "4_E3_Ecoli", "25_D1_Ecoli", "25_E3_Ecoli", "37_D1_Ecoli", "37_E3_Ecoli"]

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
drop_classes = ["LPC", "LPE", "LPE-O", "LPE-O/P", "PS-O/P"]
#drop_classes = []
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

sample_numers = ["A", "B", "C"]



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

outliers = ["4_E3_Ecoli_A", "25_D1_Ecoli_A", "25_E3_Ecoli_A", "37_E3_Ecoli_C"]

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

            exp_is = float(np.sum(exp_df_is))  # / exp_stic
            # Normalize by internal standard
            if exp_is == 0:
                print("No internal standard for class: ", uclass)
                exit()
            exp_abs = exp_df[class_bool].to_numpy()  # / exp_mtic

            exp_normed = exp_abs / exp_is

            exp_normed_sum.append(np.sum(exp_normed))

            cdf = df[class_bool]

            tl = cdf["Tail Lengths"]
            tu = cdf["Tail Unsaturation"]
            tld = cdf["Tail Length Diff"]
            tud = cdf["Tail Unsat Diff"]

            bool_tl = tl >= 0
            bool_tu = tu >= 0

            exp_avg_tl = np.average(tl[bool_tl], weights=exp_normed[bool_tl])
            exp_avg_tu = np.average(tu[bool_tu], weights=exp_normed[bool_tu])

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
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df[sample] = class_sums[i]

class_df.to_excel(writer, index=False, sheet_name="Classes")

# Save the results
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df[sample] = class_tls[i]

class_df.to_excel(writer, index=False, sheet_name="Tail Lengths")

# Save the results
class_df = pd.DataFrame(uclasses, columns=["Class Name"])
# Add results to the data frame
for i, sample in enumerate(samples):
    class_df[sample] = class_tus[i]

class_df.to_excel(writer, index=False, sheet_name="Tail Unsaturation")

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

writer.close()
