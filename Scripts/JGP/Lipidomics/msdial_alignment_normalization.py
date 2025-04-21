import pandas as pd
from Scripts.LipiDec import get_class, get_overall_tails
import matplotlib.pyplot as plt
import numpy as np
from Scripts.MTM.Lipidomics.mass_defect_analysis import kendrick_analysis
import matplotlib as mpl


def standardize_calibrant_areas(standard_names, standard_classes, samples, equisplash_volumes, eq_std_volume, split_data):
    """
    This function finds each spiked internal standard and standardizes the areas to a specific volume.

    Parameters:
    :param standard_names: List of standard names
    :param standard_classes: List of standard classes
    :param samples: List of sample numbers/identifiers
    :param equisplash_volumes: List of equisplash spike volumes
    :param eq_std_volume: Volume to which you will standardize. This is relatively arbitrary.
    :param split_data: Dataframe of the MS-DIAL output with the first 4 rows removed.

    Returns:
    2D array of standardized areas for each standard in each sample.
    """
    calibrated_areas = list([[0 for _ in range(len(samples))] for _ in range(len(standard_classes))])

    for index, row in split_data.iterrows():
        for std in range(len(standard_names)):
            if row.iloc[3] == standard_names[std]:
                for i in range(len(samples)):
                    calibrated_areas[std][i] = float(row.iloc[35 + i]) * (eq_std_volume / equisplash_volumes[i])
    return calibrated_areas

def get_normalized_areas_total_class_areas(standard_classes, standard_names, samples, sample_volumes, sample_std_volume, split_data, calibrated_areas):
    """
    This function normalizes the area of each sample-native lipid to the area of the spiked internal standard.
    Once these internal standard-normalized areas are calculated, they are subsequently normlized to the total
    lipid intensity in each sample.

    parameters:
    :param standard_classes: List of standard classes
    :param samples: List of sample numbers/identifiers
    :param sample_volumes: List of sample volumes
    :param sample_std_volume: Sample volume to which you will standardize. This is relatively arbitrary.
    :param split_data: Dataframe of the MS-DIAL output with the first 4 rows removed.
    :param calibrated_areas: 2D array of standardized areas for each standard in each sample.

    Returns:
    2D array of normalized areas for each lipid in each sample.
    1D array of total lipid intensities for each sample.
    1D array of lipid names.
    1D array of lipid classes.
    1D array of lipid adducts.
    1D array of lipid m/z values.
    1D array of lipid retention times.
    """

    normalized_intensities = [[] for _ in range(len(samples))]

    names = []
    classes = []
    adducts = []
    mzs = []
    rts = []

    for index, row in split_data.iterrows():
        #This also filters out low scoring lipids because the get_class method doesn't work on those I believe.
        name = row.iloc[3]
        lipid_class = get_class(name)
        for c in range(len(standard_classes)):
            if lipid_class == standard_classes[c]:
                if name not in standard_names:
                    names.append(name)
                    classes.append(lipid_class)
                    adducts.append(row.iloc[4])
                    mzs.append(row.iloc[2])
                    rts.append(row.iloc[1])
                    for i in range(len(samples)):
                        area = float(row.iloc[35 + i])
                        normalized_area = (area * (sample_std_volume / sample_volumes[i])) / \
                                          calibrated_areas[c][i]
                        normalized_intensities[i].append(normalized_area)

    total_lipid_intensities = [0 for i in range(len(samples))]
    total_lipid_normalized_intensities = [[] for i in range(len(samples))]

    for i in range(len(samples)):
        total_lipid_intensities[i] = np.sum(normalized_intensities[i])

        for j in range(len(normalized_intensities[i])):
            total_lipid_normalized_intensities[i].append(normalized_intensities[i][j] / total_lipid_intensities[i])

    return total_lipid_normalized_intensities, total_lipid_intensities, names, classes, adducts, mzs, rts

def get_percent_compositions(standard_classes, samples, exp_classes, total_lipid_normalized_intensities):
    """
    This function calculates the percent composition of each lipid class in each sample.

    Parameters:
    :param standard_classes: List of standard classes
    :param samples: List of sample numbers/identifiers
    :param exp_classes: List of lipid classes
    :param total_lipid_normalized_intensities: 2D array of normalized areas for each lipid in each sample.

    Returns:
    2D array of percent compositions for each lipid class in each sample.
    """
    percent_comps = [[0 for _ in range(len(standard_classes))] for _ in range(len(samples))]

    for i in range(len(samples)):
        for std in range(len(standard_classes)):
            percent_comp = 0
            for j in range(len(exp_classes)):
                if exp_classes[j] == standard_classes[std]:
                    percent_comp += total_lipid_normalized_intensities[i][j] * 100
            percent_comps[i][std] = percent_comp

    return percent_comps

def get_avg_tl_tu(samples, total_lipid_normalized_intensities, names):
    """
    This function calculates the average tail length and unsaturation for each sample.

    Parameters:
    :param samples: List of sample numbers/identifiers
    :param total_lipid_normalized_intensities: 2D array of normalized areas for each lipid in each sample.
    :param names: List of lipid names.

    Returns:
    1D array of weighted average tail lengths for each sample.
    """
    weighted_average_tail_length = list([0 for _ in range(len(samples))])
    weighted_average_unsaturation = list([0 for _ in range(len(samples))])

    for i in range(len(samples)):
        running_tl = 0
        running_tu = 0
        running_weights = 0
        for a in range(len(total_lipid_normalized_intensities[i])):
            tname, tl, tu = get_overall_tails(names[a])
            running_tl += total_lipid_normalized_intensities[i][a] * tl
            running_tu += total_lipid_normalized_intensities[i][a] * tu
        weighted_average_tail_length[i] = running_tl
        weighted_average_unsaturation[i] = running_tu

    return weighted_average_tail_length, weighted_average_unsaturation

def plot_kmds(mzs, rts, names, classes):
    kmds = [0 for i in range(len(mzs))]
    rts = [float(i) for i in rts]

    for i in range(len(mzs)):
        mass, kmass, kmd = kendrick_analysis(float(mzs[i]), 1)
        tname, tl, tu = get_overall_tails(names[i])
        kmds[i] = kmd
        plt.text(rts[i], kmd-0.001, names[i])
        if classes[i] == "CL":
            if tu == 0:
                plt.plot(rts[i], kmd, 'ro')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'rs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'r^')

        elif classes[i] == "PG":
            if tu == 0:
                plt.plot(rts[i], kmd, 'bo')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'bs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'b^')

        elif classes[i] == "PE":
            if tu == 0:
                plt.plot(rts[i], kmd, 'go')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'gs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'g^')


    plt.xlabel("Retention Time")
    plt.ylabel("KMD")
    plt.show()

def plot_tl_tu():
    samples = [1,2,3,4,5,6,7,8,9]
    tl_means = [33.6,33.97,34.34,35.945,35.095,38.905,34.325,41.1,37.945]
    tl_sds = [4.10,0.13,0.11,0.25,0.83,1.2,0.47,1.94,3.63]

    plt.bar(samples, tl_means, yerr=tl_sds, align='center', alpha=0.5, ecolor='black', capsize=10)
    plt.xticks(samples)
    plt.xlabel("Sample")
    plt.ylabel("Weighted Average Tail Length")
    plt.show()

    tu_means = [1.26,1.373,1.26,1.23,1.249,1.17,1.09,1.23,1.14]
    tu_sds = [0.32,0.04,0.03,0.13,0.31,0.07,0.15,0.13,0.09]

    plt.bar(samples, tu_means, yerr=tu_sds, align='center', alpha=0.5, ecolor='black', capsize=10)
    plt.xticks(samples)
    plt.xlabel("Sample")
    plt.ylabel("Weighted Average Unsaturation")
    plt.show()

# Set backend to Agg
mpl.use('WxAgg')

data = pd.read_table("Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\All_Reps\\Height_1_2024_08_23_14_32_45.txt", sep="\t")

#removes the various cruft stored by MS-DIAL..
split_data = data.iloc[4:]

#Update these for additional internal standards.
standard_names = ["CL 56:0|CL 28:0_28:0", "PG 33:1(d7)|PG 15:0_18:1(d7)", "PE 33:1(d7)|PE 15:0_18:1(d7)"]
standard_classes = ["CL", "PG", "PE"]

#These are the sample identifiers, equisplash volumes, and sample volumes.
#It is important that these all have the same length and order (i.e. the first sample is the first volume, etc.)
#The order and number must match the order and number of samples searched in MS-DIAL.

replicate = 2
samples = [1,4,7,2,5,8,3,6,9]
equisplash_volumes = [25,25,25,25,25,25,25,25,25]
sample_volumes = [50,50,50,50,50,50,50,50,50]

#this is arbitrary, but to be able to compare we need to standardize to a volume
#all this does is scale the final result.
standardization_volume = 50
sample_standardization_volume = 100


calibrated_areas = standardize_calibrant_areas(standard_names, standard_classes, samples, equisplash_volumes, standardization_volume, split_data)

output_df = pd.DataFrame()

total_lipid_normalized_intensities, total_lipid_intensities, names, classes, adducts, mzs, rts = (
    get_normalized_areas_total_class_areas(standard_classes, standard_names, samples, sample_volumes, sample_standardization_volume, split_data, calibrated_areas))

#percent_compositions = get_percent_compositions(standard_classes, samples, classes, total_lipid_normalized_intensities)

#avg_tl, avg_tu = get_avg_tl_tu(samples, total_lipid_normalized_intensities, names)

plot_kmds(mzs, rts, names, classes)



output_df["Lipid Name"] = names
output_df["Lipid Class"] = classes
output_df["Adduct"] = adducts
output_df["m/z"] = mzs
output_df["RT"] = rts


for i in range(len(samples)):
    output_df[str(samples[i])] = total_lipid_normalized_intensities[i]

output_df.to_excel('D:\\Johnny\\Lipidomics\\Collaborator_Samples\\Normalized_Abundances_Rep' + str(replicate) + '.xlsx', index=False)


for i in range(len(samples)):
    plt.figure(figsize=(7,7))
    plt.pie(percent_compositions[i][:], labels=standard_classes, autopct='%1.1f%%')
    plt.axis('equal')
    plt.title("Rep: " + str(replicate) +" Sample: " + str(samples[i]) + " Composition")
    plt.savefig("D:\\Johnny\\Lipidomics\\Collaborator_Samples\\Rep" + str(replicate) + "_Sample_" + str(samples[i]) + "_Composition.png")
















