import os
import matplotlib as mpl

def ensure_correct_counts(df, sample_num):

    lipids_to_remove = []
    current_lipid = "lipid"
    current_sample_idx = 0
    current_lipid_count = 0
    for i, row in df.iterrows():
        if current_lipid_count == 0:
            current_lipid = row["Molecule"]
            current_lipid_count = 0
            current_sample_idx = 0
        if row["Molecule"] == current_lipid:
            current_lipid_count += 1
            current_sample_idx += 1
            if current_sample_idx == len(sample_num):
                current_sample_idx = 0
            if current_lipid_count == len(sample_num):
                current_lipid_count = 0
        else:
            if current_lipid_count != len(sample_num):
                lipids_to_remove.append(current_lipid)
                print("Removing " + current_lipid + " from the dataset.")
                current_lipid = row["Molecule"]
                current_lipid_count = 1
                current_sample_idx = 1

    #now collect the indices of the rows with lipids to remove
    rows_to_remove = []
    for i, row in df.iterrows():
        if row["Molecule"] in lipids_to_remove:
            rows_to_remove.append(i)

    return df.drop(rows_to_remove)

def get_standardized_calibrant_areas(df, standard_names, equ_vols, equ_std_vol, sample_nums):
    '''
    This function takes Skyline quantification results as a pandas data frame,
    and returns the area of each of the spiked internal standards in each replicate.
    '''
    std_areas = [np.zeros(len(standard_names)) for i in range(len(equ_vols))]

    rows_to_remove = []
    current_sample_idx = 0
    for i, row in df.iterrows():
        name = row["Molecule"]
        if name in standard_names:
            print(row["Replicate Name"], sample_nums[current_sample_idx])
            rows_to_remove.append(i)
            std_idx = standard_names.index(name)
            std_areas[current_sample_idx][std_idx] = row["Area"] * (equ_std_vol / equ_vols[current_sample_idx])
            if np.isnan(row["Area"]):
                std_areas[current_sample_idx][std_idx] = 0
            current_sample_idx += 1
            if current_sample_idx == len(equ_vols):
                current_sample_idx = 0

    # unique_sample_nums = np.unique(sample_nums)
    # for i in unique_sample_nums:
    #     #Plot the areas of the standards in each replicate
    #     for j in range(len(standard_names)):
    #         areas = [std_areas[k][j] for k in range(len(equ_vols)) if sample_nums[k] == i]
    #         labels = np.arange(len(areas))
    #         plt.bar(labels, areas)
    #         plt.xticks(labels, labels)
    #         plt.title(standard_names[j] + " in Sample " + str(i))
    #         plt.xlabel("Replicate")
    #         plt.ylabel("Area")s
    #         plt.show()


    return std_areas, rows_to_remove

def normalize_lipid_areas(df, samp_vols, samp_std_vol, std_areas=None, std_names=None, std_classes=None, normalize_to_is=True):
    '''
    This function takes Skyline quant results as a pandas dataframe and standardizes the areas
    of each identified lipid. If standard areas are not provided, then normalization to equisplash
    will not be performed, and only sample volume and total identified area normalization will be performed.

    This function returns the input pandas dataframe with the normalized areas.
    '''

    current_sample_idx = 0

    for i, row in df.iterrows():
        name = row["Molecule"]
        if std_names is not None and name in std_names:
            continue
        else:
            lipid_class = get_class(name)
            if std_classes is not None:
                if lipid_class in std_classes:
                    class_idx = std_classes.index(lipid_class)
                    if std_areas[current_sample_idx][class_idx] == 0:
                        df.loc[i,"Area"] = 0
                    else:
                        if np.isnan(row["Area"]):
                            df.loc[i,"Area"] = 0
                        if normalize_to_is:
                            df.loc[i,"Area"] = (row["Area"] * (samp_std_vol / samp_vols[current_sample_idx])) / std_areas[current_sample_idx][class_idx]
                        else:
                            df.loc[i, "Area"] = (row["Area"] * (samp_std_vol / samp_vols[current_sample_idx])) / 1
                    current_sample_idx += 1
                    if current_sample_idx == len(samp_vols):
                        current_sample_idx = 0
            else:
                row["Area"] = row["Area"] * (samp_std_vol / samp_vols[current_sample_idx])
                current_sample_idx += 1
                if current_sample_idx == len(samp_vols):
                    current_sample_idx = 0
                if current_sample_idx == len(samp_vols):
                    current_sample_idx = 0
    return df

def write_aggregated_output(df, sample_nums, normalize_to_total=True):
    '''
    This function takes a pandas dataframe of Skyline quantification results and writes an aggregated
    output file with the mean and standard deviation of each lipid across replicates.
    '''
    total_areas = np.zeros(len(sample_nums))
    names = []
    file_names = ["" for i in range(len(sample_nums))]
    lipid_areas = [[] for i in range(len(sample_nums))]

    current_sample_idx = 0
    first_round = 1
    for i, row in df.iterrows():
        if first_round == 1:
            file_names[current_sample_idx] = row["Replicate Name"]
        if np.isnan(row["Area"]):
            row["Area"] = 0
        total_areas[current_sample_idx] += row["Area"]
        current_sample_idx += 1
        if current_sample_idx == len(sample_nums):
            current_sample_idx = 0
            first_round = 0

    current_sample_idx = 0
    for i, row in df.iterrows():
        if current_sample_idx == 0:
            names.append(row["Molecule"])
        if normalize_to_total:
            lipid_areas[current_sample_idx].append(row["Area"]/total_areas[current_sample_idx])
        else:
            lipid_areas[current_sample_idx].append(row["Area"])
        current_sample_idx += 1
        if current_sample_idx == len(sample_nums):
            current_sample_idx = 0

    output_df = pd.DataFrame()
    output_df["Lipid Name"] = names
    unique_nums = []
    for i in range(max(sample_nums) + 1):
        indices = [j for j, x in enumerate(sample_nums) if x == i]
        if len(indices) > 0:
            unique_nums.append(i)
            for k in range(len(indices)):
                column_name = "S" + str(i) + "_R" + str(k + 1)
                areas = np.array(lipid_areas[indices[k]])
                output_df[column_name] = areas
            means = []
            stds = []

            for l in range(len(names)):
                means.append(np.mean([lipid_areas[indices[k]][l] for k in range(len(indices))]))
                stds.append(np.std([lipid_areas[indices[k]][l] for k in range(len(indices))]))
            output_df["S" + str(i) + "_average"] = means
            output_df["S" + str(i) + "_StdDev"] = stds

    output_df.to_excel('Z:\\Group Share\\JGP\\CS MS1 Runs\\mzML\\normalized_abundances.xlsx',index=False)








# Set backend to Agg
mpl.use('WxAgg')

file = "Z:\\Group Share\\JGP\\CS MS1 Runs\\mzML\\CS_MS1_quant_results.csv"
os.chdir(os.path.dirname(file))

sky_df = pd.read_csv(file)



# Processing this collaborator data will be a bit difficult since each replicate was processed slightly differently.
#These numbers must be in the same order as they appear in the .csv file
sample_vols = [1.25,1.25,1.25,
               50,50,
               10,10,10,
               30,30,30,
               50,50,50,
               92,92,92,
               30,30,30,
               85,85,85,
               30,30,30,
               30,30,30,
               100,100,100,
               50,50,50,
               30,30,30,
               50,50,50,
               50,50,50,
               30,30,30,
               40,40,40,
               100,100,100,
               50,50,50,
               50,50,50,
               30,30,30,
               30,30,30,
               30,30,30,
               50,50,50,
               100,100,100,
               70,70,70]


equi_vols = [1.25,1.25,1.25,
             25,25,
             10,10,10,
             15,15,15,
             25,25,25,
             40,40,40,
             15,15,15,
             45,45,45,
             15,15,15,
             15,15,15,
             50,50,50,
             25,25,25,
             15,15,15,
             25,25,25,
             25,25,25,
             15,15,15,
             40,40,40,
             50,50,50,
             25,25,25,
             25,25,25,
             15,15,15,
             15,15,15,
             15,15,15,
             25,25,25,
             50,50,50,
             40,40,40]


sample_num = [2,2,2,
              2,2,
              8,8,8,
              2,2,2,
              9,9,9,
              7,7,7,
              1,1,1,
              4,4,4,
              5,5,5,
              3,3,3,
              3,3,3,
              8,8,8,
              4,4,4,
              7,7,7,
              5,5,5,
              8,8,8,
              5,5,5,
              1,1,1,
              3,3,3,
              6,6,6,
              9,9,9,
              7,7,7,
              6,6,6,
              4,4,4,
              6,6,6,
              9,9,9]



sky_df = ensure_correct_counts(sky_df, sample_num)

#Update these for additional internal standards.
standard_names = ["CL 56:0|CL 28:0_28:0", "PG 33:1(d7)|PG 15:0_18:1(d7)", "PE 33:1(d7)|PE 15:0_18:1(d7)"]
standard_classes = ["CL", "PG", "PE"]

#first find the standards in each replicate
std_areas, rows_to_remove = get_standardized_calibrant_areas(sky_df, standard_names, equi_vols, 50, sample_num)

sky_df = sky_df.drop(rows_to_remove)

sky_df = normalize_lipid_areas(sky_df, sample_vols, 50, std_areas, standard_names,
                               standard_classes, normalize_to_is=False)

write_aggregated_output(sky_df, sample_num, normalize_to_total=True)

sky_df.to_csv("normalized_skyline_quant_results.csv", index=False)

