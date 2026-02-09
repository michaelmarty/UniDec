import os
import pandas as pd
from sec_plotting import folder_handler, save_param, style_settings
from matplotlib import pyplot as plt
import numpy as np
import gc


def spec_plotter(array_dict , colours=None, xlim=None, ylim=None, title="Normalized Native Mass Spectra", grouping_index=None, subplots=False):
    """
    This function compiles a bunch of plot settings and whatnot to generate native and CD mass spec figures. Differs from sec plotting because it does
    :param array_dict: the dictionary containing the pandas dataframe(s) containing your mass spec data
    :param colours: a list of colours to be applied to each line/group
    :param xlim: optional tuple of limits to use for the x-axis
    :param ylim: optional tuple of limits to use for the y-axis should contain at most 2 floats
    :param title: optional title of the figure; str defaults to normalized native mass spectra
    :param grouping_index: List; defaults to None; User defined grouping factor
    :param subplots: Bool; defaults to False; Forces a stacked plot of data
    :return: None
    """
    xlab = "Mass (kDa)"
    ylab = "Relative Intensity (%)"
    x = "Mass"
    y = "Intensity"
    style_settings()
    counter = 0
    if subplots:
        if grouping_index:
            length = len(grouping_index)
        else:
            length = len(array_dict)//3
        subplotter(array_dict, length, ylab, colours=colours, title="Normalized Native Mass Spectra", grouping_index=grouping_index)
    else:
        for spectra in array_dict:
            plt.plot(array_dict[spectra][x], array_dict[spectra][y], label=spectra, c=colours[counter])
            counter += 1
        plt.legend(loc=2, frameon=False, labelspacing=0.1)
        plt.title(title)
        plt.ylabel(ylab)

    if xlim:
        plt.xlim(xlim)
    else:
        plt.xlim(0)
    if ylim:
        plt.ylim(ylim)
    else:
        plt.ylim((0,100.01))
    plt.xlabel(xlab)
    plt.yticks(ticks=[0, 50, 100], labels=['0', '%', '100'])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    save_bool = input("Do you want to save the figure? (y/n)")
    if save_bool == "y":
        plt.savefig(save_param(), transparent=True)
    plt.show()

def subplotter(array_dict, length, ylab, grouping_index=None, colours=None, title=None):
    fig, axs = plt.subplots(length, sharey=True, figsize=(32, 20))
    fig.subplots_adjust(hspace=.7)
    midd = length // 2
    axs[midd].set_ylabel(ylab)
    fig.align_ylabels(axs)
    counter = 0
    x = "Mass"
    y = "Intensity"
    if grouping_index:
        keys = list(array_dict.keys())
        for index in grouping_index:
            axs[counter].plot(array_dict[keys[index]][x], array_dict[keys[index]][y], label=keys[index], c=colours[index])
            axs[counter].legend(loc=1).get_frame().set_alpha(0)
            counter += 1
    else:
        rep = 1
        for chroma in array_dict.keys():
            axs[counter].plot(array_dict[chroma][x], array_dict[chroma][y], label=chroma, c=colours[counter])
            axs[counter].legend(loc=2).get_frame().set_alpha(0)
            if counter > 0:
                axs[counter].sharex(axs[0])
            if rep == 3:
                axs[counter].spines['top'].set_visible(False)
                axs[counter].spines['right'].set_visible(False)
                counter += 1
                rep = 1
            else:
                rep += 1

    axs[0].set_title(title)

def spec_normalized(array_dict):
    """
    Normalizes the mass spectra so that the relative intensities can be compared across different species
    :param array_dict: the dictionary containing the pandas dataframe(s) containing your mass spec data
    :return: the dictionary containing the altered pandas dataframe containing normalized mass spectra
    """
    for spectra in array_dict:
        array_dict[spectra]["Intensity"] = 100*array_dict[spectra]["Intensity"] / np.max(array_dict[spectra]["Intensity"])
        array_dict[spectra]["Mass"] = array_dict[spectra]["Mass"] / 1000
    return array_dict

def main():
    """
    Use this function to create plots and figures. Documentation for other functions in this file can be found
    in the README.md file.
    :return:
    """

    print("All Done thank you!")

if __name__ == "__main__":
    main()
    gc.collect()
