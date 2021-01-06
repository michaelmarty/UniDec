UniDec: Universal Deconvolution of Mass and Ion Mobility Spectra 
=================================================================

UniDec is a Bayesian deconvolution program for deconvolution of mass spectra and ion mobility-mass spectra.

It was originally published in: [M. T. Marty, A. J. Baldwin, E. G. Marklund, G. K. A. Hochberg, J. L. P. Benesch, C. V. Robinson, Anal. Chem. 2015, 87, 4370-4376.](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00140)

Detailed descriptions of the algorithm are provided in the paper and subsequent papers. Please cite us if you use UniDec in your research.

Please contact mtmarty@email.arizona.edu for questions, suggestions, or with any bugs.

## Installation

UniDec may be downloaded from [https://github.com/michaelmarty/UniDec/releases](https://github.com/michaelmarty/UniDec/releases).

This compiled version is compatible with 64-bit Windows. It a portable binary, so it does not need a conventional installation.
Just unzip the folder, put it somewhere convenient, and click the GUI_UniDec.exe file in the folder to launch.

## Tutorial

You can watch a video tutorial on how to use UniDec here: [https://www.youtube.com/watch?v=e33JxgY6CJY](https://www.youtube.com/watch?v=e33JxgY6CJY).

## Licensing

UniDec is distributed under a completely open source license. Our hope is that this allows UniDec to be
more widely used. If you are interested in including UniDec in another academic or commercial software distribution, 
you are welcome to email mtmarty@email.arizona.edu for more information. 

UniDec source code and compiled binaries are released under a modified BSD License as described below. Note, we ask
that you cite us in any publications. Quantitative citation metrics will help grant applications to support future development.

By downloading UniDec, you are agreeing to the UniDec and any third party licensing agreements and all terms therein. 

### UniDec License:

Copyright (c) 2016, University of Oxford
              2017-2019, University of Arizona
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions, and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.
4. Any publications that result from use of the software should cite Marty et al. Anal. Chem. 2015. DOI: 10.1021/acs.analchem.5b00140. If UniDec is redistributed or incorporated into other software, it must be clearly indicated to the end user that UniDec is being used, and the request to cite Marty et al. Anal. Chem. 2015. DOI: 10.1021/acs.analchem.5b00140 must be passed on to the end user.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### Third Party Licenses

Waters DLLs are distributed under the Waters MassLynxSDK EULA license. The Thermo RawFileReader DLLs are distributed under the RawFileReader.doc license. By downloading and using these DLLs, you are accepting those licenses and these terms. Thermo and Waters are indemnified, defended, and held harmless from any claims, including attorneys’ fees, related to the distribution or use of this software. Redistribution of Waters and Thermo libraries is restricted as described in those licenses. Info on other licenses are provided below.

## UniDec Compatible File Types

UniDec is built to open .txt files using [numpy.loadtxt](http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.loadtxt.html). 

For MS data, it opens a two-column either a tab or space delimited list of m/z and intensity values.

For IM-MS, it will open a three-column tab or space delimited list of m/z, arrival time (or bin), and intensity values. Sparse matrices are fine for IM-MS. 

It is compatible with a text header at the beginning of the file. It will skip lines until it reaches the start of the data.

For Water's .raw files, UniDec is bundled with converters (CDCReader.exe) to 
convert the data to .txt. It will compress the retention time dimension into a single spectrum. 
A single file can be opened directly, or multiple files can be converted using 
Tools > Simple Batch Process Raw to Txt. For a fancier conversion such as extracting specific functions or scans, 
try Tools > Raw to Txt Conversion Wizard. Note: rawreader.exe has been replaced with the MassLynxSDK 4.5 Python library. Water's converters will need MassLynxRaw.dll and/or cdt.dll in the same directory as the converter executables (the unidec_bin folder or the top directory). You can find these at: [https://interface.waters.com/masslynx/developers-area/sdks/](https://interface.waters.com/masslynx/developers-area/sdks/). 

Thermo .raw files can be read as you would a text file on Windows thanks to [multiplierz](https://github.com/BlaisProteomics/multiplierz). You will need [MSFileReader](https://thermo.flexnetoperations.com/control/thmo/download?element=6306677) installed. Please cite them (http://onlinelibrary.wiley.com/doi/10.1002/pmic.201700091/abstract). It will compress all scans together unless parsed with MetaUniDec. 

Note: multiplierz is currently not compatible for Python 3. However, you can make it work by fixing the print and import commands. I can provide a modified version on request.

Finally, many vendor formats can be converted mzML using [Proteowizard](http://proteowizard.sourceforge.net/). UniDec will open mzML file as if they are a text file, and this format should be cross platform.
We utilize [pymzML](http://pymzml.github.io/intro.html#general-information) for this. Please [cite them](https://www.ncbi.nlm.nih.gov/pubmed/22302572).

## MetaUniDec File Types and Importing

With MetaUniDec, everything is stored in a single HDF5 files. 
The HDF5 Import Wizard allows you to import a range of different file types directly into a single HDF5 file.
Thermo RAW and mzML files are supported fully, which means that either the scan or time range can be specified.
Text and Waters RAW files are supported for file import. Text files must be a single m/z spectrum.
Waters RAW files will have all scans summed into a single m/z spectrum upon import. 
The File>Waters Conversion Wizard tool allows specific scans to be converted into text files for importing.

In addition to the Import Wizard, there are several Manual File options, which will allow you to create blank HDF5 
(New File) and load data into it (Add File). Note, Add Data will sum all scans together, and Waters data is not supported.
You can select multiple files at once here. 
You can also just copy the data from XCalibur or MassLynx and then use Add Data From Clipboard. 

There are a few automated tools to parse single chromatograms directly into HDF5 files if you have all the data chromatograms 
with predictable scans or times. You can batch process multiple files at once. 
Only mzML and Thermo RAW files are supported for automated chromatogram import.

## Installing the Source Code

Most users will likely just want to run the compiled version. For those advanced users who have experience with Python,
we have provided the source code for the GUI and API.

### Python

UniDec is currently compatible only with Python 3. There are several Python libraries that UniDec will depend on. 

matplotlib
numpy
scipy
wxpython
natsort
pymzml
networkx
h5py
pypubsub
tornado
pythonnet
multiplierz (Windows only, for Agilent imports)

All of these can be installed from the command line with (for example):
    
    pip install natsort

Note: I would highly recommend setting up 64-bit Python as the default. MS data works fine with 32-bit, but IM-MS data is prone to crash the memory. If you are getting memory errors, the first thing to try is to upgrade the bit level to 64.

### The UniDec Binaries

As described below, the Python code presented here relies on one critical binary, UniDec.exe. The binary should be in the /unidec_bin directory. 

If you are interested in building the binary or modifying the source code, the code and Visual Studio project files
are in unidec_src/UniDec. It is currently configured for Visual Studio Community 2015 with HDF5 1.10.1 and Intel Parallel Studio 19.
It can be easily compiled with other compilers but will show a significant performance loss without the Intel Compilers.
 
If you are interested in binaries for Mac and Linux, I would recommend building them yourself using the scripts in the unidec_src/UniDec directory. UniDec compiles easily and works fine with these operating systems, but I don't have the time to support them.

## UniDec Documentation

Documentation is for the Python engine and GUI and can be found at http://michaelmarty.github.io/UniDecDocumentation/.

My goal is that this documentation will allow you to utilize the power of the UniDec python engine for scripting data analysis routines and performing custom analysis. Also, it will allow you to add new modules to the UniDec GUI.

I'm still working on documenting some of the windows and extensions (including MetaUniDec and C code), but the core features should be here.

## UniDec Architecture


UniDec is bilingual. The core of the algorithm is written in C and compiled as a binary.
It can be run independently as a command line program fed by a configuration file.

The Python engine and GUI serve as a very extensive wrapper for the C core. 

The engine (unidec.py) can be operated by independently of the GUI. This allows scripting of UniDec analysis for more complex and high-throughput analysis than is possible with the GUI.
The engine contains three major subclasses, a config, data, and peaks object.

The GUI is organized with a Model-Presenter-View architecture.
The main App is the presenter (GUniDec.py).
The presenter contains a model (the UniDec engine at unidec.py) and a view (mainwindow.py). 
The presenter coordinates transfer between the GUI and the engine.

MetaUniDec has a similar structure with a presenter (mudpres.py), engine (mudeng.py), and view
 (mudview.py). However, unlike conventional UniDec, MetaUniDec includes a number of additional features to process 
 and analyze data. It relies less heavily on the Python API.


## Getting Started


Here is some sample code for how to use the engine. 

    import unidec
    
    file_name="test.txt"
    folder="C:\\data"
    
    eng=unidec.UniDec()
    
    eng.open_file(file_name, folder)
    
    eng.process_data()
    eng.run_unidec(silent=True)
    eng.pick_peaks()

In reading the documentation, it is perhaps best to start with the unidec.UniDec class.
The main GUI class is GUniDec.UniDecApp.

## Change Log

v.4.4.0

**Added Scroll Bar to Controls!** People have been asking for this for a while, and I could never figure out how to get it to work. Finally, I managed to find the answer and added it with some buttons on the bottom to help expand and collapse key parts.

In a significant technical change, I switched all doubles in the C code to floats. For those who are interested, this sacrifices a little precision in calculations for improvements in speed and file size.

**Added Fast Profile and Fast Centroid** options for UniChrom2 and MetaUniDec. These will help speed up deconvolutions and limit file sizes by not adding the massgrid and mzgrid to each spectrum. The animate features will not work, and other things might not be available, but the basic settings should be consistent.

**Major speed improvements to large mzML data files** by using gzip to compress the data prior to opening it. Thanks to the pymzML team for this.

**Added DScore to MetaUniDec and UniChrom**. This will now show the DScore for each peak. It calculates the DScore for each spectrum but then takes the average across all spectra weighted by the intensity of that peak in each spectrum. 
Also, added Filter by Score to Experimental Menu.

**Added Sliding Window to UniChrom**. You can now specify the width of the window (in minutes) and the offset between the start of the windows (in # of scans). The offset needs to be an integer greater than or equal to 1. Setting the window to 1 will start a window on each scan. Setting a window of 0 and an offset of 1 will give every scan separately, without averaging any of them together.

**Hidden feature: Write data from plot to text file.** Clicking Ctrl+u on most plots will now give you a dialog to save the underlying data as a text file. Usually, this data was written somewhere behind the scenes, but this will give an easy way to export it.

Added DNA to Biopolymer Calculator.

Added a Mass Defect Comparison Tool to the Mass Defect Window.

Added normalization option to FFT window, which really improves the results.

Added "Subtract Constant" option to the UniDec baseline subtraction. The number in the box will specify the constant relative to the max intensity. For example, a value of 0.01 will subtract each intensity by 1% of the maximum.

UniChrom2 manual selection now puts data out as text files in the UniDec Files and Figures folder rather than TestSpectra.

Bug fixes to data import functions.

Fixed bugs with figure saving and isolating peaks in UniChrom.

v.4.3.0

**Added UniChrom2: UniDec for LC/MS data**. UniChrom is built on top of the MetaUniDec engine and uses HDF5 files and many of the same core tools. The primary additions are the ability to parse chromatography data into HDF5, visualize and interact with LC/MS data, and manually select and quickly deconvolve parts of the chromatogram.

**Added Data Reduction in MetaUniDec and UniChrom**. This now mirrors the behavior of UniDec for processing data to remove a fixed percentage of the lowest intensity data.

**Added Native Combining of Waters Data**, which should dramatically speed up averaging Waters chromatograms by using the native MassLynx libraries.

**Added Thermo RawFileReader libraries** to avoid having to install MSFileReader and to make opening native Thermo data faster and more robust.

The mass defect windows can now make horizontal lines from input mass values, not just input mass defect value. It will automatically calculate the mass defect from the input mass.

Major refactoring of the code to support UniChrom. For example, switched from local to absolute paths for most files in the engine.

Added Estimated Area extraction to DataCollector and MetaUniDec. Here, it uses the peak height and FWHM to estimate the area based on the peak shape.

Added experimental DoubleDec (Double Deconvolution) feature to the main UniDec window.

Fixed bug with spaces in oligomer names when importing ofiles.

Switched default native charge range to +/- 1000 from +/- 100 to work better with denatured data.

Other bug fixes and minor changes.

v.4.2.2

**Added Open Recent Files** menu item.

**Added Copy All Full** to copy paste a wide range of parameters from the peaks in the peak list. Let me know if there are additional parameters you would like added.

**Major upgrades to file mzML and Thermo file imports**. The Thermo files should now open and average faster. Both types should now allow large files to open without crashing the memory.

Added ability to use a semicolon with the X values in the Data Collector to grab the sum of multiple states. Also, added a feature to Plot X Ranges to help to visualize the integration range.

Added time midpoint to the mzML parsing.

Added the version number to the config file for future use or reference.

Added experimental DoubleDec feature to Data Collector. 

Added experimental report generator for the Native MS Guided Structural Biology Center.

Fixed bug with fits being plotted when process data is clicked. Fixed a number of additional bugs. Hopefully didn't create too many more...

v.4.2.1

**Added a Full button** to reset the m/z range to the full value.

**Added menu item to "Load Prior State for Current File"**. This will load the past deconvolution results, so you don't have to click deconvolve again if you've already run the data. This now has taken the Ctrl + L keyboard shortcut, which was previously Load State. Load and Save State have been renames to Load and Save Zip State File for clarity.

Improved reliability of auto peak width function for sparse data. 

Fixed issues with plotting and deconvolution with interpolation between sparse data points.

Fixed issue with deleting spectra in MetaUniDec.

Fixed issue where markers were clipped on the edges of plots. This introduces some issues with single data points not being clipped when spectra are zoomed in. I'm working to fix this for future versions.

Updates to add fitting to the Mass Defect Extractor.

v.4.2.0

**Added Smart Transform**. From the beginning, UniDec has asked users to decide how it should perform the nonlinear transform from m/z to mass (Integration or Interpolation). Now, the Smart Transform decides for the user based on the density of the data. Where the data is sparse relative to the mass sampling, it will use interpolation to fill in the gaps. Where data is dense relative to the mass sampling, it will use integration to ensure that all the data in-between is accounted for. Smart Transform is located under the Additional Deconvolution Parameters tab and has now been made the default. In some preliminary testing, the Smart Transform showed more robust peak areas that were less sensitive to the mass sampling. The peak heights may change some, but the peak areas should be a lot more reliable (remember Ctrl+I to integrate). Note: you may not notice much difference, but using the Smart Transform should avoid some glitches that are hard to spot. Let me know how it goes!

**Major improvements in 2D plotting speed** due to behind-the-scenes changes and optimzation of the plotting code.

**Sparse data support**. Improvements to the algorithm to support sparse data, including just peak centroids. This fundamentally changes the underlying algorithm in minor ways. You may notice minor differences but likely will not. Let me know if it causes any problems.

Better handling of click and zoom at the edges of the plot. Now, it will default to the edge if you go off the plot. 

Added drag and drop for loading the state from a zip file.

Upgraded pymzml version for improved mzML compatibility and made mzML import more error tolerant.

Fixed issue with discrete plot and data that was non-uniformly sampled.

Fix to bug in calculating the error from the weighted standard deviation of charge state masses. 

Fixed bug in Load and Save State with non-text files.

Fixed IM-MS plotting bugs.

v.4.1.2

**Added button in UniDec to turn off data normalization.** Note: the beta values will be automatically scaled to match this. 

Renamed the parameter ZScore in the UniScore calculation to CSScore. Added R squared to the UniScore calculation.

Added several experimental subtract and divide features. Tweaks to linear regression experimental feature. 

Bug fixes and compatibility updates.

v.4.1.1

Added right click feature to color peaks by their scores.

Improvements to Data Collector module to fit protein oligomerization KDs and improve customization. 

Added experimental Linear Regression features for analysis of repeating mass units.

Added experimental FPOP feature to print and copy average degree of oxidation.

Various syntax optimizations to address deprecation warnings and library updates.

v.4.1.0

Introduced **UniScore**. This will automatically score the quality of your peaks and the overall deconvolution. Various experimental features were also added with this to visualize and filter the scores. Some updates to the scoring algorithm from Version 4.0.2 to improve it.

A good heuristic for UniScore and DScore values is:

80-100: A (Excellent)

60-80: B (Good)

40-60: C (Fair)

20-40: D (Poor)

0-20: F (Almost certainly noise)

Added **Calculator for getting masses from Protein/Peptide or RNA sequences**. This can be launched stand alone from the Experimental Menu or by right clicking on the Oligomer and Mass Tools Oligomer Builder box.  

Fixed bug to now allow protein-protein oligomerization KD values to be measured with the Data Collector. Find out more on the [UniDec Wiki Page](https://github.com/michaelmarty/UniDec/wiki/Fitting-KDs-with-Data-Collector).

v.4.0.2

Switched UniDec to load processed data if it already exists.

**Added Copy All Basic** to the peak panel, which copies the peak mass, height, and area to a format you can paste into excel. I can add additional copy modes easily now, so let me know what you would like to be able to copy out. 

Added control of spectra color map on MetaUniDec.

Added Waterfall Plots on Experimental Menu in MetaUniDec.

Many Bug Fixes: Sped up quick control responses by only doing auto peak width when needed. Fixed memory leak with HDF5 files. Fixed bugs slowing down file imports. Fixed bug with PDF Report Generator. Fix bug with isolate/ignore/repopulate/reorder in MetaUniDec.

Experimental features for peak scoring. More to come in future releases. 

v.4.0.1

Added experimental support for **Agilent .D** files. Please test it and let me know how it works. You may need to run UniDec as an administrator first to get it to register all the DLLs correctly.

Improved the peak panel by switching to white text when the background gets dark, adding commas to the masses (hopefully it will be easier to read), and **adding new right click features to label masses**. Also fixed text colors on MetaUniDec spectrum panel.

Moved the Bin Every and Background Subtraction parameters from the main Data Processing tab to the Advanced Data Processing Parameters tab. **Added a Backround Subtraction check box** to turn on background subtraction with a curved background set to 100. 

Added peak centroids, FWHM, and error between charge states to the _peakparams.dat output.

Added a check box to ignore zero values in UltraMeta and fixed the bar charts on UltraMeta to have the same colors as the main window. 

Added the ability in UltraMeta to right click on a file and open it in MetaUniDec.

Added the ability to reorder spectra in MetaUniDec by editing the indexes. 

Added new extraction choices for area of peaks above a threshold in UltraMeta and MetaUniDec.

UltraMeta will now be tolerant of different files having different x values when creating error bars. It will only plot and average consensus data points. 

Mass Defect plots will now switch to kDa from Da if appropriate.

When adding data to HDF5 files, having CID_n or SID_n in the file name will import n as the collision voltage. 

Fixed several bugs/issues to make everything run more smoothly.

v.4.0

Added **Quick Controls** to the main panel. This should allow you to turn on and off features quickly without
those pesky numbers. The advanced controls are still available as before. 

A new experimental feature has been added based to use a SoftMax function (controlled by the parameter beta)
to **suppress deconvolution artifacts**. The higher beta is, the more the algorithm with seek a single charge state assignment for each data point. Setting beta to zero will turn this off. This seems to work best when combined with other assumptions, such as peak width, point smoothing, and charge and mass smoothing. It also seems to work best when combined with background subtraction. Play around with it and let me know what you think. A minus flag applies the SoftMax to the entire data set rather than just a single column of the m/z vs charge matrix at once. 

Changed the Waters MS import from rawreader.exe to the MassLynxSDK 4.5. You may need to download new MassLynxRaw.dll files from Waters. I am working with Waters to get approval to distribute the DLL files bundled with UniDec, but it isn't final yet. Because I used the Waters API, I am not releasing the source code for the Waters Importer until I have that agreement in place. I will do this as soon as I can. In the mean time, Waters import features will work on the compiled binary version but not from the Python source code. You can keep using v.3.2 in the mean time. Stay tuned...

Added **Ctrl+C to copy out images from plots**. You should be able to paste these into other applications.

Added **Example Data**, which can be quickly loaded from the File menu. You can also add data to this by dropping your own files in the Example Data folder. It works in the same was as the custom presets.

Added **Data Reduction** data processing feature in UniDec for removing noise from large data sets. Basically, you set a percentage of the data you want to remove. UniDec then finds what intensity threshold is required to remove that much and takes out all data below that threshold. 

Added **UniChrom** for quick viewing of chromatogram TICs. You can open mzML or Thermo .Raw files directly. Waters .Raw files can be opened by dragging and dropping it in the main window.

Added more functionality for UltraMeta.

Added drag and drop for _conf.dat files to UniDec to more easily import settings.

Adjustments to the algorithm to improve speed and reliability. Switched build to Visual Studio 2019 and Intel Parallel Studio 19. 

Changes to the Gaussian blur functions for charge and mass, which are activated by negative flags for mass and charge smooth width.

Removed cross validation feature because I'm pretty sure no one was using it.

v.3.2

Added a logo to the start screen.

Update to the latest mzML specification.

Build update to Python 3.7 and latest libraries.

Bug fixes.

v.3.1.1

Expanded isotope mode to either output monoisotopic masses or average masses. 

Added an experimental feature to plot average charge state for peaks in UniDec. 

Bug fixes.

v.3.1.0

Added new parameter for Point Width Smooth. 

Changed the right control panels to streamline the key parameters. 

Updated the default parameters to make them more universal.

Added custom presets. You can drop any _conf.dat file in the Presets folder and it will add a custom preset for it. This folder can be organized into subfolders to collect your presets. A few custom presets from the Marty Lab are included. Feel free to send me yours if you would like them uploaded to the public distribution.

Note: Several background changes to the algorithm (MS-only) have allowed the use of more general settings, such as a peak width of 0 (v.2.7.3) and a point smooth width of 1, which adds an additional smooth of +/- 1 data point. We hope that these defaults and changes to the layout will allow new users to more easily access the software. Power users will still have access to these parameters in the advanced settings windows. 

v.3.0.2

Modified the FFT window tool on MetaUniDec to be the dual FFT of each spectum rather than the sum. Indiviual spectra FFT windows are still available by right clicking the spectra in the list on the left. 

Added button for Negative Ion Mode in Additional Deconvolution Parameters. This simply switches the adduct mass (typically a proton) to negative when clicked. 

Bug fixes to MetaUniDec and UltraMeta.

v.3.0.1

Added averagine isotope distribution calculator for peaks in Experimental. Thanks to Jim Prell for developing the Fourier isotope distribution calculation function.

Fixed bug with FWHM calculation.

Added in peak centroid for intensity within FWHM when FWHM is calculated.

v.3.0.0

**Updated everything to Python 3.6.**

Improvements to Mass Defect window, including new extractor.

v.2.7.3

Added experimental features for charge state smoothing when Charge or Mass Smooth Width is greater than 1. 
Added experimental feature to allow for zero peak width.
Added Kendrick Mass Defect shortcut with Ctrl+K. Switched previous Ctrl+K shortcut (plot peaks) to Ctrl+J.

v.2.7.2

Adding registerInterfaces command in when multiplierz fails for Thermo file direct reads.
For this to work, UniDec may need to be run as an administrator. 

v.2.7.1

Fixed Bug in Waters Import Wizard and UltraMeta Mass Defect plots.

v.2.7.0

Added 2D plots in Mass Defect tools for MetaUniDec.

v.2.6.8

Fixed bug with matching spectra in Oligomer and Mass Tools.

v.2.6.7

*Space+Middle Click* on any line plot will now automatically add peak annotation. Simply Space+Middle Click again to turn it off.
For those on a laptop, Alt+Left Click will also toggle between labelling and not. 

v.2.6.6

*Shift+Middle Click* on any plot will now spawn a window to specify the y range manually.
*Alt+Middle Click* on any plot will now spawn a window to specify the x range manually.
Fixed legend in UltraMeta. 

v.2.6.5

Added *Ctrl+ Middle Click* on any plot window to bring up a dialog to change the rcParams from matplotlib (https://matplotlib.org/users/customizing.html).
Added -peaks flag in UniDec.exe to get peaks for each individual spectrum in an HDF5 file.

v. 2.6.4

A few bug fixes.
Added file name as metadata from Import Wizard.
Added automatic monomer dimer assign for the Wysocki Group.
Added annotated m/z and mass animations for MetaUniDec.
Added ability to save figures automatically from animation.

v. 2.6.3

Added new Marker Selector to Peak List right click menu.
Fixed color selection bug.

v. 2.6.2

Cool new plots on the Mass Defect Tools.

v. 2.6.1

Added limits to the number of plots MetaUniDec will plot before it will plot only a representative number.
This significantly speeds up plots for very large data sets.

Sped up raw and mzML file imports.

v. 2.6.0

HDF5 Import Wizard Built for MetaUniDec. Waters can be easily imported directly through the HDF5 Import Wizard.

Updated Help Docs in MetaUniDec

Moved UltraMeta to Analysis Menu from Experimental.

Added UltraMeta and HDF5 Import Wizard to Launcher.

Added threading and fixed bugs in UltraMeta.

v. 2.5.0

Fits for FFT Window. Also, Moved FFT Window in UniDec to Analysis Menu from Experimental. 

v. 2.4.0

Update to Ultrameta Data Collector to specify the location of the peaks explicitly.

v. 2.3.0

Compatibility upgrades to support the most recent version of many of the major packages. Upgrade to wxpython 4.0, matplotlib 2.1, and scipy 1.0.

v. 2.2.0

Added **help documentation to MetaUniDec**. These should be useful for learning how MetaUniDec and its tools work.

Added importing multiple chromatograms by range of times or scans. 
This allows you to compile certain timepoints/scans from multiple files into one HDF5 file easily.

Added errors for peaks. Three types of error have been added: FWHM (both), duplicates (MetaUniDec), and mean (UniDec).

Added bar graphs to visualize the different parameters for exponential decay, line, or sigmoid fitting in UltraMeta.

Added a compare tool to the FFT Window. Once clicked, you drag boxes around the regions you want to compare, and then hit the compare button again to plot the regions.

Adjusted how baselines are calculated. UniDec has parameters added under Experimental->Additional Parameters which allow you to adjust the baseline calculation.

Added a **repack tool** to fix a known problem with HDF5 files. If you were to use Data Processing or UniDec Parameters that made the HDF5 very large, changing the parameters to values that would usually make the HDF5 small would not result in a shrinking of the HDF5 file. The Repack Directory tool will recursively repack all HDF5 files in a directory to their current size.


v. 2.1.1

Added fitting to exponential decay, line, or sigmoid in MetaUniDec and UltraMeta.

v. 2.1.0

**Added support for native opening of Thermo Raw files** using multiplierz (https://github.com/BlaisProteomics/multiplierz).
You should now be able to open RAW files directly as you would a text file and to parse the chromatograms with MetaUniDec. 
Other file types from other vendors should be possible if you send me an example file to test.

v. 2.0.0

**MetaUniDec** added with a number of tools for batch processing. This is a whole new program dedicated to expand the UniDec core to high-throughput dataset processing and visualization. It is twice as fast as conventional batch processing with UniDec and combines the dataset into a single HDF5 file.

**Launcher** added to manage the various tools. 

**UltraMeta** data collector added to analyze sets of datasets produced by MetaUniDec.

Numerous behinds the scenes changes to simplify code and make compatible with MetaUniDec. Other bug fixes and improvements.

v. 1.4.0

Added common masses table and features to Oligomer and Mass Tools.

v. 1.3.1 

Experimental peak fitting in mass defect window.

v. 1.3.0

Added automatic baseline feature in crude form. 

Updated preset configs and added a Nanodisc preset.

Removed zip save file creation from batch processing. Let me know if you would like this brought back, but it seems to be just creating unnecessary files.

v. 1.2.5

Fixed glitch with speedy flag. Removed it entirely from Python code. Linflag will now import and overwrite the speedy option in UniDec.exe if present in config file.

v. 1.2.4

Updated builds for MacOS and Linux. Updated commandline printout.

v. 1.2.3

Improvements to automatic peak width detection for nonlinear data.

v. 1.2.2

Tweaks to the resolution settings to adjust for different monitor size. Slightly increased default figure size.

Minor fixes and improvements.

v. 1.2.1

**Added Undo/Redo buttons and keyboard shortcuts to Experimental Menu.**

New extraction modes on DataCollector (50% and 10% Thresholded Center of Mass).

Fixed major bug in Load State, potential bug in PDF report, and other random bugs.

v. 1.2.0

Largely behind the scenes changes. Added HDF5 integration for C and Python but haven't committed to this as a save file yet.

Merged C code for UniDecIM.exe into UniDec.exe, so there is a single executable for both 1D and 2D.

A number of other bug fixes, updates, and subtle improvements.

v. 1.1.0

**Added Linear and Power Law calibrations for T-Wave IM-MS.** These are highly untested so proceed with caution. Please let me know how they work.

Linear: Reduced CCS = Calibration Parameter 1 * Reduced Drift Time + Calibration Parameter 2

Power Law: Reduced CCS = Calibration Parameter 1 * (Reduced Drift Time ^ Calibration Parameter 2)

(For reference, the previous log calibration was and is)
Log: Reduced CCS = Exp(Calibration Parameter 1 * log(Reduced Drift Time) + Calibration Parameter 2)

Dr. Tim Allison updated CDCReader.exe, which converts Waters Raw IMMS files into txt to perform more accurately with small bin sizes and with ending m/z cutoffs.

**This version requires an updated binary of UniDecIM.exe to utilize the new calibrations. I am still working on the IP for this, so contact me if you need this binary.**

v. 1.0.13

Added **log and square root intensity scales** under the advanced menu. Note, these don't save with the file but carry over from the session.
Added DPI control for figure outputs in save figure dialog.
Improved GUI to better save user input.
Added an advanced menu item to **open the current save directory** in the file explore of the OS.

v. 1.0.12

Added drag and drop in the data collector utility. Can drop either text files to add them to the list or a JSON save file to load it.

v. 1.0.11

Thanks to Tim for bug fixes on CDCReader.exe (the binary that converts Water's IMMS files to text files).

Updates to UniDecIM.exe to change how the twaveflag is determined and handled. Will now accept multiple possible twaveflags to allow for alternative calibration strategies in the future. Contact me with your function of interest.

Updated contact email for MTM.

v. 1.0.10

Added preset defaults for the **Exactive EMR** under "High-resolution Native". The previous high-resolution preset is now "Isotopic Resolution". 

v. 1.0.9

Small change in C code so that the peaks are now multiples of the mass bin size. Will need the updated [binary](http://unidec.chem.ox.ac.uk/).

Cleaned up single click mass calculations on the raw data plot so that a zoom out is not registered as a click. Control + click will register a click without zooming out.

v. 1.0.8

**Fixed Waters IM-MS conversion errors!** 
Thanks to Tim Allison for teaming up to track down the source of these problems.
It should now work for any IM-MS file, but it will require 64-bit in everything.

v. 1.0.7

Added a **mass tolerance on matching**. Now, matches that fail to come within a certain tolerance will be ignored. 

v. 1.0.6

Added an experimental window for grid deconvolution.
This is close to the original implementation of UniDec back when it was a Nanodisc-only script in Mathematica (see Marty, et. al. JASMS 2014).
It has come a long way, but it is still pretty crude.

v. 1.0.5

Added directory name to file names window.
Added smarter labeling of plots in data collector.
Added ability to specify Variable 1 units in data collector and add this as an axis label.
Added support for viridis and other new color maps in matplotlib 1.5.0.

v. 1.0.4

Fixed bugs in parallel processing for Data Collector Utility and Import Wizard.

v. 1.0.3

Added Waters .raw file open and .raw to .txt batch conversion to drag and drop.

v. 1.0.2

Bug fix on KD fit parallel processing.

v. 1.0.1

Added **drag and drop** on the main window. Drag a single file into the plotting area to open it. Drag multiple files to run them in batch mode.

v. 1.0

Total rewrite of the program to separate the GUI and engine. 
The engine allows scripting of UniDec in Python without the GUI.
Added documentation and cleaned up code with major refactoring and simplification.
Most changes are in the back end and will hopefully be invisible to people using the GUI.

* New Features:
    * Added **automatic peak width** determination with Ctrl + W or Tools > Automatic 
    * Added **"Display Mass Differences"** to right click menu in the peak panel. This will display the mass differences between each peak and the selected peak. Very useful (in my humble opinion).
    * Left clicking on two peaks in the m/z spectrum will **solve for the mass from the two m/z values**.
    * To allow left clicking on zoomed region, holding ctrl while clicking will prevent rescaling the axes.
    * Right click on the m/z spectrum will determine the **max and min m/z value from the current zoom on the plot**.
    * Zeroing the zoomed region is now a double right click on the m/z spectrum (useful for eliminating noise spikes).
    * **Middle click** on any plot now opens a **save figure dialog** to allow you to save that specific figure.
    * Sped up autocorrelation, added it directly to the Analysis menu, and added it to the Data Collector.
    * Added "Save Figures as .pdf" as a save figure shortcut.
    * Moved "Get Spctrum from Clipboard" to File menu and added shortcut at Ctrl + G.
    * Moved Integrate/Interpolate option for converting from m/z to mass from the Advanced menu to the Additional Filter/Restraints control window.  
    
* From the last update:
    * Report Center of Mass: Shows the center of mass for the zoomed region in the mass distribution (Plot 2)
    * Plot by Charge: Peak are now charge states rather than mass species. Plot charge state distribution. Plots each charge state as a separate distribution in Plot 4.
    * Plot Charge Offsets: Replot the Mass vs. Charge plot as a Mass vs. Charge Offset. Note: may be slow for large arrays.
    * Auto Match Tools: Basically the same as clicking Oligomer and Mass Tools > Match to Mixed Oligomers > OK. Requires the oligomers to already be defined in the Oligomer and Mass Tools window.
    * Kendrick Mass Analysis: Tools for high-mass mass defect analysis. The paper on this is in press.
    * 2D Grid Extraction: Extract intensity values for predefined mass values in a 1D or 2D grid. 

* Experimental Features (Unpublished, untested, and unfinished):
    * Calibration Window: Allows a polynomial calibration to be applied to the m/z data. It will not modify the data file or save the result.
    * FFT Window: Views double FFT for windowed regions of the spectrum. I know this is weird, but this is a sneak peak into something we are working on.
    * Color Plots: Plot specific regions of the spectrum in the color of the peak. This is common for native MS papers, but it can be a mess when there are a lot of overlapping peaks.
    * Get Errors: Working on a automated error determination...


## Additional License Info

h5py: http://docs.h5py.org/en/stable/licenses.html

wxpython: https://wxpython.org/pages/license/index.html

numpy: https://www.numpy.org/license.html

scipy: https://www.scipy.org/scipylib/license.html

matplotlib: https://matplotlib.org/users/license.html

natsort: https://github.com/SethMMorton/natsort/blob/master/LICENSE

pymzml: https://github.com/pymzml/pymzML/blob/master/LICENSE.txt

networkx: https://networkx.github.io/documentation/networkx-1.10/reference/legal.html

multiplierz: https://github.com/BlaisProteomics/multiplierz/blob/master/LICENSE

pypubsub: https://pypubsub.readthedocs.io/en/v4.0.3/about.html#license

License-Zoombox and Zoomspan:

Copyright (c) 2010, Duke University and the United States Department of
Veterans Affairs. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of Duke University nor the United States Department
     of Veterans Affairs may be used to endorse or promote products derived
     from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS
IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

For portions of this code, copyright and license information differs from
the above. In these cases, copyright and/or license information is inline.






