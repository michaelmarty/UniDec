from matplotlib import pyplot as plt, rcParams
import numpy as np
import pandas as pd
import os
import myqueue as mq
import gc


def style_settings():
    """
    Style settings for the marty lab called in later functions
    :return: None
    """
    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['lines.linewidth'] = 1
    rcParams['errorbar.capsize'] = 3
    rcParams['patch.force_edgecolor'] = True
    rcParams['patch.facecolor'] = 'b'
    rcParams['lines.markersize'] = 7
    rcParams['font.size'] = 18
    rcParams['font.weight'] = 500
    rcParams['font.sans-serif'] = "Arial"
    rcParams['axes.titleweight'] = 500
    rcParams['axes.titlesize'] = 22
    rcParams['axes.linewidth'] = 3
    rcParams['xtick.major.width'] = 1
    rcParams['ytick.major.width'] = 1



def save_param():
    """
    User input to make a file path for saving/writing files
    :return: master_pth: str; a path to be passed into plt.savefig()
    """
    file_name = input("Enter file name: ")
    # file_pth = input("Enter file path (include '\\'): ")
    # file_name = "CDMS_masterplot"
    file_pth = "Z:\\Group Share\\Estevan\\Aiden\\Replicates\\Discs - done\\Plots\\"
    master_pth = file_pth + file_name
    return master_pth


def plot_chrom(array_dict: dict,normalized=False,xlim=None, ylim=None, grouping_index=None,subplots=False,colours = None, title = "SEC Chromatograms"):
    """
    Much like spec_plotter() this function takes chromatograms and plots them all on one graph
    :param array_dict: the dictionary containing each array (chromatogram) to be plotted Keys == names
    :param normalized: bool; if the data has been normalized this updates the ticks to match
    :param xlim: optional x-axis limits for smaller windows
    :param ylim: optional y-axis limits for shorter windows
    :param grouping_index: optional grouping parameter to plot only select spectra from the dictionary
    :param colours: lst; a list of colors associated with each line must match in length the number of names
    :param title: title for the graph; defaults to 'SEC Chromatograms'
    :return: None (plots)
    """
    style_settings()
    xlab = "Retention Volume (mL)"
    if normalized:
        ylab = "Normalized Absorbance"
    else:
        plt.ylim(0, np.max(array_dict))
        ylab = "Absorbance (mAU)"
    if subplots:
        if grouping_index:
            length = len(grouping_index)
        else:
            length = len(array_dict)
        subplotter(array_dict, length, ylab, grouping_index, colours, title)
    else:
        no_subplots(array_dict, grouping_index, colours)
        plt.ylabel(ylab)
        plt.title(title)
        plt.tight_layout()
    if xlim:
        plt.xlim(xlim)
    else:
        plt.xlim(0, 40)
    if ylim:
        plt.ylim(ylim)
    plt.yticks(ticks=[0, .5, 1], labels=['0', '0.5', '1'])
    plt.ylim(-0.01, 1)
    plt.xlabel(xlab)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    save_bool = input("Do you want to save the figure? (y/n)")
    if save_bool == "y":
        plt.savefig(save_param())
    plt.show()

def subplotter(array_dict, length, ylab, grouping_index=None, colours=None, title=None):
    fig, axs = plt.subplots(length, sharey=True, figsize=(16, 10))
    fig.subplots_adjust(hspace=.7)
    midd = length // 2
    axs[midd].set_ylabel(ylab)
    fig.align_ylabels(axs)
    counter = 0
    x = "ml"
    y = "mAU"
    if grouping_index:
        keys = list(array_dict.keys())
        for index in grouping_index:
            axs[counter].plot(array_dict[keys[index]][x], array_dict[keys[index]][y], label=keys[index], c=colours[index])
            axs[counter].legend(loc=1).get_frame().set_alpha(0)
            counter += 1
    else:
        for chroma in array_dict.keys():
            axs[counter].plot(array_dict[chroma][x], array_dict[chroma][y], label=chroma, c=colours[counter])
            axs[counter].legend(loc=1).get_frame().set_alpha(0)
            if counter > 0:
                axs[counter].sharex(axs[0])
            counter += 1
    axs[0].set_title(title)

def no_subplots(array_dict, grouping_index=None, colours = None):
    counter = 0
    x = "ml"
    y = "mAU"
    if grouping_index:
        keys = list(array_dict.keys())
        for index in grouping_index:
            plt.plot(array_dict[keys[index]][x], array_dict[keys[index]][y], label=keys[index], c=colours[index])
    else:
        for chroma in array_dict.keys():
            plt.plot(array_dict[chroma][x], array_dict[chroma][y] ,label=chroma, c=colours[counter])
            counter += 1
    plt.legend(loc="center right", frameon=False, labelspacing=0.05)


def normalize_sec(array_dict):
    """
    Normalizes each absorbance height for each chromatogram passed into it
    :param array_dict: sec data with named label: Pandas array
    :return: array_dict: dict: normalized absorbance array
    """
    for chroma in array_dict.keys():
        maxx = np.max(array_dict[chroma]["mAU"])
        array_dict[chroma]["mAU"] = array_dict[chroma]["mAU"] / maxx
    return array_dict


def folder_handler(folder_name: str, names: list = None, formatted="csv"):
    """
    This function reads in data saved in either csv or txt files and saves into a workbook of sorts
    make sure that the length of names provided matches the number of files or you'll get an error
    :param folder_name: str; name of the folder holding the data
    :param names: lst; the names you want associated with the data must be in folder order (check in PyCharm)
    :param formatted: str; this is to denote whether you're reading in a .csv or .txt file defaults to csv
    *note* .txt files are ALWAYS treated as mass spectra
    :return: array_dict: dict{str:pd.dataframe}; a dictionary containing the dataframes mapped to the names provided
    """
    file_queue = mq.Queue()
    counter = 0
    array_dict = {}
    try:
        any(file_queue.enqueue(file) for file in os.listdir(folder_name))
        while not file_queue.is_empty:
            file_path = os.path.join(folder_name, str(file_queue.dequeue()))
            if formatted == "csv":
                temp = pd.read_csv(file_path, header=2, usecols=[0, 1])
            if formatted == "txt":
                temp = pd.read_csv(file_path, names=["Mass", "Intensity"], delimiter=" ")
            array_dict[names[counter]] = temp
            del temp
            counter += 1
    except FileNotFoundError:
        print('Please Ensure the Folder Exists')
    except UnicodeDecodeError:
        print('Please make sure every file is in the correct format')
    except UnboundLocalError:
        print("Please use a supported file format")
    return array_dict


def main():
    """
    This function handles all of the plotting change only and uses the above functions to do so
    For your own plots you will need to update....
    names to contain a list of the names corresponding to your folder
    color to contain a list of the colors you want
    folder_name in the first function call to be the folder your data is in
    grouping_index in whichever function you need to plot specific lines
    title in EVERY plot function to match accurately the data you're plotting
    :return: Plots out the wazoo
    """

    print("All Done! Thanks")
    gc.collect()



if __name__ == "__main__":
    main()
