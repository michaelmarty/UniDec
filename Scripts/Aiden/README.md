---
output: html_document
toc: True
---

# Getting Started

Hello this file should serve as a guide to using the functions included in `spec_plotter.py` and `sec_plotting.py`
I'll try and keep it updated as new plots/features get added


As of 11.16 this has not been updated to fit changed made since the last update

This README was updated 1.27.26 by Aiden Wolf to match updates in functions and file names

## User Facing Functions:

### `folder_handler():`

#### Parameters: 

`folder_name` : the folder where the data files are hosted; need file path if not in working directory

`names` : a list of the names of the data you're reading in

`formatted` : a str input that denotes the type of file you're reading in; defaults to csv; accepts: csv and txt

***Note:** txt files will automatically be treated as mass spectra (as of 11.7.25)*

#### Returns: 

`array_dict`: this is a dictionary that holds every data frame you imported with the function and is stored such that 
`dict = {name: dataframe}` 

This function aims to read through a folder provided by the user as `folder_name=` to bring each data file into python. 

#### It errors if:
- The number of names does not match the number of files
- The file format does not match the `formatted=` argument

#### To fix errors:

- Ensure the number of names you've input matches the number in the directory
- Make sure you've saved the file correctly (this is mainly true for chromatograms)

### `normalize_sec()` and `spec_normalized()`

#### Parameters:

`array_dict`: This can be passed straight from `folder_handler` and its generally good practice **not** to call it before you're finished plotting all of your normal plots

#### Returns:

`array_dict`: super helpful I know. Essentially it returns the same dictionary of arrays that it receives but with a normalized y-axis so that everything is relative.

These functions both essentially do the same thing, however for mass spec the mass gets divided by 1,000 and the intensity gets normalized to a percentage.
In the SEC function the absorbance gets normalized to a proportion

This function family usually does not error. 


## Functions that handle plotting

I made these originally as separate function simply because they started as very divorced projects, but now I've refused to put them together simply because it's easier to leave defaults in place. I might combine them at some point into one big plotting script... until then...


### `plot_chrom()`

#### Parameters:

`array_dict`: the dictionary containing each array (chromatogram) to be plotted Keys == names

`normalized`: bool; if the data has been normalized this updates the ticks to match

`xlim`: can be a tuple `()` or a list `[]` of one or two limits that you want to use as a cut off for the axes

`ylim`: same as xlim

`colours`: a list `[]` of colours. Should match, in length, the number of chromatograms you're plotting

`title`: Optional title to be given to the plot. It defaults to "SEC Chromatograms"


#### Errors:

Most errors in this function come from the other functions not working... but one big one can come if your colours list is smaller than the number of files you have

### `spec_plotter()`

#### Parameters:

`array_dict`: the dictionary containing each array (chromatogram) to be plotted Keys == names

`colours`: a list `[]` of colours. Should match, in length, the number of spectra you're plotting

`xlim`: can be a tuple `()` or a list `[]` of one or two limits that you want to use as a cut off for the axes

`ylim`: same as xlim

`title`: Optional title to be given to the plot. Defaults to Normalized Native Mass Spectra

`grouping_index`: List; an optional grouping term that can be used to plot only specific spectra; defaults to `None`

`subplots`: Bool; tells the plotter to make stacked plots by feeding the data into the `subplotter()` function. Subplots are grouping index compatible; defaults to `False`

#### Errors:

Most errors in this function come from the other functions not working... but one big one can come if your colours list is smaller than the number of files you have


## Non-User Facing Functions:

### `style_settings()`

This function has no parameters and returns nothing it is simply a container for all the style settings for publishable figures
Any alterations should be limited. This is called *always* and is applied to every graphical output. You should ONLY change this if Dr. Marty says you need something in a different style

### `save_param()`

**This technically a user facing function as users interact with it.**

Essentially this function serves as a method to save figures with names specified by the user without having to directly alter a function

The function asks for a figure name and file path. The file path can be empty in which case it writes into the working directory

### `subplotter()`

This function makes the plot fit into a stacked subplot. The MS version currently is configured for replicate data, but I'll add an option for non replicates soon... (1.27.26)

_Note that there's two separate versions of this function_

**USERS SHOULD ALTER THE FIGURE SIZE DEPENDING ON NUMBER OF SPECTRA**

### `no_subplots`

Handles plotting a chromatogram on one plot. Will be depreciated soon