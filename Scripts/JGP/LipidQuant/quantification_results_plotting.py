import os
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd

from Scripts.LipiDec.AddFeaturestoDF import *

is_names = ["CL 56:0|CL 28:0_28:0",
            "PI 15:0|18:1(d7)",
            "PG 33:1(d7)|PG 15:0_18:1(d7)",
            "SM 18:1|18:1(d9)",
            "PC 15:0|18:1(d7)",
            "PE 33:1(d7)|PE 15:0_18:1(d7)",
            "PS 15:0|18:1(d7)"]

def parse_file_id(id):
    split = id.split("_")
    sample = int(split[0][1:])
    br = int(split[1][2:])
    tr = int(split[2][2:])
    return sample, br, tr

def plot_class_piecharts(spectra_files, normalized_quant_results):
    os.chdir("Z:\\Group Share\\JGP\\Lipidomics\\Skyline_testing\\LipidQuant_piecharts")
    sample_contents = {}

    for file in spectra_files:
        class_areas = {}
        for i in range(len(normalized_quant_results[file.fileid])):
            lipid = normalized_quant_results["Lipid Name"][i]
            lipid_class = get_class(lipid)
            if lipid_class not in class_areas:
                class_areas[lipid_class] = normalized_quant_results[file.fileid][i]
            else:
                class_areas[lipid_class] += normalized_quant_results[file.fileid][i]

        if file.sample_num in sample_contents:
            if file.biorep in sample_contents[file.sample_num]:
                sample_contents[file.sample_num][file.biorep][file.techrep] = class_areas
            else:
                sample_contents[file.sample_num][file.biorep] = {}
                sample_contents[file.sample_num][file.biorep][file.techrep] = class_areas
        else:
            sample_contents[file.sample_num] = {}
            sample_contents[file.sample_num][file.biorep] = {}
            sample_contents[file.sample_num][file.biorep][file.techrep] = class_areas

    #Now need to calculate the percentage of each class within each sample number,
    #biological replicate, and technical replicate

    sample_class_areas = {}
    for sample_num in sample_contents:
        sample_class_areas[sample_num] = {}
        for biorep in sample_contents[sample_num]:
            class_fractions = {}
            for techrep in sample_contents[sample_num][biorep]:
                total_content = np.sum(list(sample_contents[sample_num][biorep][techrep].values()))
                for lipid_class in sample_contents[sample_num][biorep][techrep]:
                    if lipid_class in class_fractions:
                        class_fractions[lipid_class].append(sample_contents[sample_num][biorep][techrep][lipid_class]/total_content)
                    else:
                        class_fractions[lipid_class] = [sample_contents[sample_num][biorep][techrep][lipid_class]/total_content]

            #This gets the average class area within each biological replicate
            average_class_areas = {}
            for lipid_class in class_fractions:
                average_content = np.average(list(class_fractions[lipid_class]))
                average_class_areas[lipid_class] = average_content

            sample_class_areas[sample_num][biorep] = average_class_areas

            labels = list(average_class_areas.keys())
            fractions = list(average_class_areas.values())

            plt.pie(fractions, labels=labels, autopct='%1.1f%%')
            plt.title("Sample:" + str(sample_num) + " BR:" + str(biorep))
            plt.axis('equal')
            plt.savefig("Sample" + str(sample_num) + "_BR" + str(biorep)+".png", dpi=300)
            plt.show()
    return

def plot_indlipid_barcharts(spectra_files, normalized_quant_results):
    samples_bioreps_techreps = {}
    for file in spectra_files:
        if file.sample_num in samples_bioreps_techreps:
            if file.biorep in samples_bioreps_techreps[file.sample_num]:
                samples_bioreps_techreps[file.sample_num][file.biorep].append(file.techrep)
            else:
                samples_bioreps_techreps[file.sample_num][file.biorep] = [file.techrep]
        else:
            samples_bioreps_techreps[file.sample_num] = {}
            samples_bioreps_techreps[file.sample_num][file.biorep] = [file.techrep]

    for i in range(len(normalized_quant_results["Lipid Name"])):
        lipid_name = normalized_quant_results["Lipid Name"][i]
        #Plot the area in each techrep of each biorep in each sample
        for sample_num in samples_bioreps_techreps:
            br_areas = {}
            for biorep in samples_bioreps_techreps[sample_num]:
                tr_areas = [0.0 for tr in range(len(samples_bioreps_techreps[sample_num][biorep]))]
                for j in range(len(samples_bioreps_techreps[sample_num][biorep])):
                    tr = samples_bioreps_techreps[sample_num][biorep][j]
                    id  = str("S" + str(sample_num) + "_BR" + str(biorep) + "_TR" + str(tr))
                    tr_areas[j] = normalized_quant_results[id][i]

                br_areas[biorep] = tr_areas

            brs = list(br_areas.keys())
            avgs = [np.mean(vals) for vals in list(br_areas.values())]

            plt.bar(brs, avgs)
            plt.title("Sample:" + str(sample_num) + " Lipid:" + lipid_name)
            plt.show()

def compare_skyline_lipidquant(spectra_files, normalized_quant_results):
    skyline_df = load_skyline_results()
    skyline_ratio_dict = get_ratio_dict(skyline_df)
    lipidquant_ratio_dict = get_ratio_dict(normalized_quant_results)

    replicate_ids = []
    for column in skyline_df.columns.tolist():
        if "Lipid" in column or "Standard" in column:
            continue
        else:
            replicate_ids.append(column)


    for id in replicate_ids:
        ratio_pairs = []


        max_ratio  = 0
        classes_ratios = defaultdict(list)
        for lipid in normalized_quant_results["Lipid Name"]:
            lipid_class = get_class(lipid)
            if lipid in skyline_ratio_dict and lipid in lipidquant_ratio_dict:
                skyline_ratio = skyline_ratio_dict[lipid][id]
                if skyline_ratio > max_ratio:
                    max_ratio = skyline_ratio
                lipidquant_ratio = lipidquant_ratio_dict[lipid][id]
                if lipidquant_ratio > max_ratio:
                    max_ratio = lipidquant_ratio
                ratio_pairs.append([skyline_ratio, lipidquant_ratio])
                classes_ratios[lipid_class].append((skyline_ratio, lipidquant_ratio))

        #Now make a scatterplot of the ratio pairs

        # Assign a color to each class automatically
        for lipid_class, ratios in classes_ratios.items():
            x_vals = [r[0] for r in ratios]
            y_vals = [r[1] for r in ratios]
            plt.scatter(x_vals, y_vals, label=lipid_class)
        plt.title(id)
        plt.legend(title="Lipid Class", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.xlabel("Skyline Ratio")
        plt.ylabel("Lipid Quant Ratio")
        plt.xlim(left=0, right=max_ratio*0.1)
        plt.ylim(bottom=0, top=max_ratio*0.1)
        plt.show()
    return

def load_skyline_results():
    skyline_df = pd.read_csv(r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\RawSkyline_MS1_Quantification.csv")

    rep_fileid_dict = {8:"S1_BR1",
                       25: "S1_BR2",
                       13: "S1_BR3",
                       26: "S3_BR1",
                       21: "S3_BR2",
                       24: "S3_BR3",
                       17: "S2_BR1",
                       15: "S2_BR2",
                       9: "S2_BR3"}

    #The purpose of this is to load these results into a format similar
    #to what is spit out by LipidQuant
    quantified_lipids = {}

    lipid_names = []
    lipid_classes = []
    lipid_int_stds = []
    lipid_name_idx_dict = {}

    added_lipids = 0

    for index, row in skyline_df.iterrows():
        lipid_name = row["Molecule"]
        lipid_class = row["Molecule List Name"]
        is_internal_std = lipid_name in is_names

        if lipid_name not in lipid_names:
            lipid_names.append(lipid_name)
            lipid_classes.append(lipid_class)
            lipid_int_stds.append(is_internal_std)
            lipid_name_idx_dict[lipid_name] = added_lipids
            added_lipids += 1

    replicate_areas = {}
    for id in rep_fileid_dict:
        for tr in range(1,4):
            file_id = rep_fileid_dict[id] + "_TR" + str(tr)
            replicate_areas[file_id] =  [None] * len(lipid_names)

    for index, row in skyline_df.iterrows():
        lipid_name = row["Molecule"]
        idx = lipid_name_idx_dict[lipid_name]
        replicate = row["Replicate Name"]
        area = float(row["Area"])


        base_id = rep_fileid_dict[replicate]
        for tr in range(1,4):
            possible_id = base_id + "_TR" + str(tr)
            if replicate_areas[possible_id][idx] is not None:
                continue
            else:
                replicate_areas[possible_id][idx] = area
                break

    for id in replicate_areas:
        replicate_areas[id] = [0.0 if x is None else x for x in replicate_areas[id]]

    skyline_lipidquant_df = pd.DataFrame()
    skyline_lipidquant_df["Lipid Name"] = lipid_names
    skyline_lipidquant_df["Lipid Class"] = lipid_classes
    skyline_lipidquant_df["Standard"] = lipid_int_stds
    for id in replicate_areas:
        skyline_lipidquant_df[id] = np.array(replicate_areas[id])

    return skyline_lipidquant_df

def get_ratio_dict(quant_results):
    replicate_ids = []
    for column in quant_results.columns.tolist():
        if "Lipid" in column or "Standard" in column:
            continue
        else:
            replicate_ids.append(column)

    # These will eventually be removed from the quantification dataframe, but not yet.
    standard_class_indices = {}
    for i in range(len(quant_results["Standard"])):
        if bool(quant_results["Standard"][i]):
            standard_class_indices[quant_results["Lipid Class"][i]] = i

    # First level of this is class, second is sample number
    # third level is BR which contains a list of non-zero areas
    sample_br_standard_areas = {}
    for lipid_class in standard_class_indices:
        if lipid_class not in sample_br_standard_areas:
            sample_br_standard_areas[lipid_class] = {}
        for rep_id in replicate_ids:
            sample, br, tr = parse_file_id(rep_id)
            std_area_in_file = quant_results[rep_id][standard_class_indices[lipid_class]]
            if std_area_in_file > 0:
                if sample in sample_br_standard_areas[lipid_class]:
                    if br in sample_br_standard_areas[lipid_class][sample]:
                        sample_br_standard_areas[lipid_class][sample][br].append(std_area_in_file)
                    else:
                        sample_br_standard_areas[lipid_class][sample][br] = [std_area_in_file]
                else:
                    sample_br_standard_areas[lipid_class][sample] = {}
                    sample_br_standard_areas[lipid_class][sample][br] = [std_area_in_file]


    sample_br_std_averages = {}
    for lipid_class in sample_br_standard_areas:
        if lipid_class not in sample_br_std_averages:
            sample_br_std_averages[lipid_class] = {}
        for sample in sample_br_standard_areas[lipid_class]:
            sample_br_std_averages[lipid_class][sample] = {}
            for br in sample_br_standard_areas[lipid_class][sample]:
                sample_br_std_averages[lipid_class][sample][br] = np.mean(sample_br_standard_areas[lipid_class][sample][br])

    ratio_dict = {}
    for i,row in quant_results.iterrows():
        lipid_name = row["Lipid Name"]
        lipid_class = get_class(lipid_name)
        ratio_dict[lipid_name] = {}
        for id in replicate_ids:
            area_in_file = row[id]
            std_area_in_file = quant_results[id][standard_class_indices[lipid_class]]
            if std_area_in_file > 0:
                ratio = area_in_file / std_area_in_file
                ratio_dict[lipid_name][id] = area_in_file / std_area_in_file
            else:
                sample, br, tr = parse_file_id(id)
                std_area_in_file = sample_br_std_averages[lipid_class][sample][br]
                ratio = area_in_file/std_area_in_file
                ratio_dict[lipid_name][id] = area_in_file / std_area_in_file

    return ratio_dict






