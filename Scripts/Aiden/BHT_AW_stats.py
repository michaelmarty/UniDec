import myqueue as mq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def folder_handler(folder_name: str, und=8):
    """
    Function that iterates through a folder of .txt (or .dat) files from UniDecChromCD and returns an array of dictionaries.
    Assumes that the file name has the injection number after 8 sets of underscores... will fix eventually
    :param folder_name:
    :return: array_dict: dict{key: pd.DataFrame}; a dictionary of arrays pooled from all of the files in the folder
    """
    file_queue = mq.Queue()
    counter = 0
    array_dict = {}
    try:
        any(file_queue.enqueue(file) for file in os.listdir(folder_name))
        while not file_queue.is_empty:
            file_path = os.path.join(folder_name, str(file_queue.dequeue()))
            number = file_path.split("_")[und]
            temp = pd.read_csv(file_path, delimiter=" ", header=0, usecols=["Mass", "Area", "MassFWHM"])
            key = int(number)
            array_dict[key] = temp
            del temp
            counter += 1
    except FileNotFoundError:
        print('Please Ensure the Folder Exists')
    except UnicodeDecodeError:
        print('Please make sure every file is in the correct format')
    except UnboundLocalError:
        print("Please use a supported file format")
    except IndexError:
        print("fix me")
    return array_dict

def comparison(replicates: list, var="all", vers='sliding'):
    """
    Does comparisons of each ion statistic and writes results for each into a file/group of files depending on
    :param replicates:
    :param var: str;
    :param vers: str; the version of comparison test required. Two options: sliding and standard. Sliding compares local variation whereas standard compares each value to the highest number of injections.
    :return:
    """
    if vers == 'sliding':
        pass
    elif vers == 'standard':
        pass
    pass


def tolerant_index(array_dict, ref_mass=0, threshold=0):
    lower = ref_mass - threshold
    upper = ref_mass + threshold
    first_form = {}
    key_list = sorted(list(array_dict.keys()))
    for key in key_list:
        for i in range(len(array_dict[key]["Mass"])):
            if lower < array_dict[key]['Mass'][i] < upper:
                first_form[key] = i
                break
    del key_list, upper, lower, ref_mass
    return first_form


def merge_replicates(replicates: list, n=3, select=None, length = 0):
    """
    Averages all of the injections for each replicate in a set and returns a single dictionary of averaged arrays with standard deviation
    :param replicates: list; list of array dictionaries to be merged
    :param n: int; number of replicates
    :return: array_dict: dict{key: pd.DataFrame}; key is the injection number and pd.DataFrame is the averaged array with standard deviation for each variable
    """
    array_dict = {}
    key_list = sorted(list(replicates[0].keys()))
    if select:
        for key in key_list:
            for i in range(len(select)):
                start = select[i][key]
                end = select[i][key] + length
                replicates[i][key] = replicates[i][key][start:end]
    for key in key_list:
        averaged_vars_frame = pd.DataFrame(
            np.mean([replicates[0][key], replicates[1][key], replicates[2][key]], axis=0),
            columns=replicates[0][1].columns)
        standard_dev_frame = pd.DataFrame(
            np.std([replicates[0][key], replicates[1][key], replicates[2][key]], axis=0, ddof=1) / np.sqrt(3),
            columns=["SEM Mass", "SEM Area", "SEM MassFWHM"])
        avg_frame = averaged_vars_frame.merge(standard_dev_frame, left_index=True, right_index=True)
        array_dict[key] = avg_frame
    return array_dict


def ppm(ref, actual):
    """
    Calculates the ppm difference between two values
    :param ref: int; reference value
    :param actual: int; current value for comparison
    :return: ppm difference; int
    """
    return ((ref-actual)/ref)*(10**6)

def converge_plot(array_dict: dict, ions=False,colors=None, ensemble=None):
    data = {}
    error = {}
    keys = list(array_dict.keys())
    nrows = len(array_dict[keys[0]]["Mass"])
    fig, ax = plt.subplots(nrows,  sharey=True, figsize=(16, 10))
    for i in range(len(array_dict[1]["Mass"])):
        for key in keys:
            if not ions:
                y = np.abs((array_dict[len(array_dict)]["Mass"][i] - array_dict[key]["Mass"][i]))
                yerr = np.abs(array_dict[len(array_dict)]["SEM Mass"][i] - array_dict[key]["SEM Mass"][i])
                data[key] = y
                error[key] = yerr

        ax[i].scatter(data.keys(), data.values())
        ax[i].plot(data.keys(), data.values())
        ax[i].errorbar(x=data.keys(), y=data.values(), yerr=error.values(),fmt='o' ,ecolor='black', elinewidth=2, capsize=5, capthick=2)
    plt.show()

def main():
    pmr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 1")
    pmr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 2")
    pmr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 3")
    replicates = [pmr1, pmr2, pmr3]
    pmavg = merge_replicates(replicates)
    converge_plot(pmavg)
    del pmr1, pmr2, pmr3, pmavg
    # Since there is an increase in inherent heterogeneity- will add a tolerance threshold for the 5 most abundant glycoforms and apply statistical testing to those.
    # For starters, will comb through until the first glycoform is found within 50 Da and the following 4 rows will be included for comparisons until Rabbit mAbs are used
    wmabr1 =  folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 1", und=9)
    wmabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 2", und=9)
    wmabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 3", und=9)
    select_r1 = tolerant_index(wmabr1, 148220, 20)
    select_r2 = tolerant_index(wmabr2, 148220, 20)
    select_r3 = tolerant_index(wmabr3, 148220, 20)
    select = [select_r1, select_r2, select_r3]
    wmabavg = merge_replicates([wmabr1, wmabr2, wmabr3], select=select, length=5)
    converge_plot(wmabavg)

if __name__ == "__main__":
    main()