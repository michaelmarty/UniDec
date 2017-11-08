UniDec: Universal Deconvolution of Mass and Ion Mobility Spectra 
=================================================================

UniDec is a Bayesian deconvolution program for deconvolution of mass spectra and ion mobility-mass spectra.

It was orignally published in: [M. T. Marty, A. J. Baldwin, E. G. Marklund, G. K. A. Hochberg, J. L. P. Benesch, C. V. Robinson, Anal. Chem. 2015, 87, 4370-4376.](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00140)

Detailed descriptions of the algorithm are provided in the paper. Please cite us if you use UniDec in your research.

Please contact mtmarty@email.arizona.edu for questions, suggestions, or with any bugs.

##Installation

UniDec may be downloaded from [https://github.com/michaelmarty/UniDec/releases](https://github.com/michaelmarty/UniDec/releases).

This compiled version is compatible with 64-bit Windows. It a portable binary, so it does not need a conventional installation.
Just unzip the folder, put it somewhere convenient, and click the GUI_UniDec.exe file in the folder to launch.

##Licensing

We have recently converted to a completely open source license. Our hope is that this allows UniDec to be
more widely used. If you are interested in including UniDec in another academic or commercial software distribution, 
you are welcome to email mtmarty@email.arizona.edu for more information. 

UniDec source code and compiled binaries are released under a modified BSD License as described below. Note, we ask
that you cite us in any publications. Quantitative citation metrics will help grant applications to support future development.

UniDec License:

Copyright (c) 2016, University of Oxford
              2017, University of Arizona
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
4. Any publications that result from use of the software must cite Marty
   et al. Anal. Chem. 2015. DOI: 10.1021/acs.analchem.5b00140.

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

## UniDec Compatible File Types

UniDec is built to open .txt files using [numpy.loadtxt](http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.loadtxt.html). 

For MS data, it opens a two-column either a tab or space delimited list of m/z and intensity values.

For IM-MS, it will open a three-column tab or space delminated list of m/z, arrival time (or bin), and intensity values. Sparse matrices are fine for IM-MS. 

Recent versions are compatible with a text header at the beginning of the file. It will skip lines until it reaches the start of the data.

For Water's .raw files, UniDec is bundled with converters (CDCReader.exe and rawreader.exe) to 
convert the data to .txt. It will compress the retention time dimension into a single spectrum. 
A single file can be opened directly, or mulitiple files can be converted using 
Tools > Simple Batch Process Raw to Txt. For a fancier conversion such as extracting specific functions or scans, 
try Tools > Raw to Txt Conversion Wizard.

Water's converters will need [MassLynxRaw.dll](https://www.waters.com/waters/support.htm?locale=en_GB&lid=134815627&cid=511442&type=DWNL) (32-bit, for MS) and/or [cdt.dll](https://www.waters.com/waters/support.htm?locale=en_GB&lid=134825112&cid=511442&type=DWNL) (64-bit, for IM-MS) in the same directory as the converter executables (the unidec_bin folder or the top directory).

Thermo .raw files can be read as you would a text file on Windows thanks to [multiplierz](https://github.com/BlaisProteomics/multiplierz). You will need [MSFileReader](https://thermo.flexnetoperations.com/control/thmo/download?element=6306677) installed. Please cite them (http://onlinelibrary.wiley.com/doi/10.1002/pmic.201700091/abstract). It will compress all scans together unless parsed with MetaUniDec. 

Finally, many vendor formats can be converted mzML using [Proteowizard](http://proteowizard.sourceforge.net/). UniDec will open mzML file as if they are a text file, and this format should be cross platform.
We utilize [pymzML](http://pymzml.github.io/intro.html#general-information) for this. Please [cite them](https://www.ncbi.nlm.nih.gov/pubmed/22302572).

## MetaUniDec File Types

With MetaUniDec, everything is stored in a single HDF5 files. 
There are a few automated tools to parse chromatograms into HDF5 files if you have all the 
data chromatograms with predictable scans or times. Otherwise, you need to load the data in. 
First, create a new HDF5 using New File. Then, use Add Data Files to add text files, Thermo RAW, 
or mzML files into the HDF5 file. You can select multiple files at once here. 
You can also just copy the data from XCalibur or MassLynx and then use Add Data From Clipboard. 
That is probably the easiest way to get started and then you can work on optimizing the workflow from there. 
You can use the Waters converter if you need to get waters data into text files, 
and proteowizard can convert a variety of vendor files to mzML.

## Installing the Source Code

Most users will likely just want to run the compiled version. For those advanced users who have experience with Python,
we have provided the source code for the GUI and API.

### Python

UniDec is currently compatible only with Python 2.7. There are several Python libraries that UniDec will depend on. 

matplotlib==1.5.1
numpy
scipy
wxpython==3.0.2
natsort
twython
pymzml
networkx
h5py
multiplierz (Windows only, for Thermo RAW imports)

With the exception of wxpython, which can be installed from the [web](http://wxpython.org/), all of these can be installed from the command line with (for example):
    
    pip install natsort

Note: I would highly recommend setting up 64-bit Python as the default. MS data works fine with 32-bit, but IM-MS data is prone to crash the memory. If you are getting memory errors, the first thing to try is to upgrade the bit level to 64.

### The UniDec Binaries

The Python code presented here relies on one critical binary, UniDec.exe. 
UniDecIM.exe has been merged into UniDec.exe and is no longer used.
The binary should be in the /unidec_bin directory. 

If you are interested in building the binary or modifying the source code, the code and Visual Studio project files
are in unidec_src/UniDec. It is currently configured for Visual Studio Community 2015 with HDF5 1.10.1 and Intel Parallel Studio 17.
It can be easily compiled with other compilers but will show a significant performance loss without the Intel Compilers.
 
If you are interested in binaries for Mac and Linux, they are also in the unidec_bin directory as unidecmac and unideclinux.
However, we do not build these regularly, so I would recommend building them yourself using the scripts in the unidec_src/UniDec directory.

##UniDec Documentation


Documentation is for the Python engine and GUI and can be found at http://michaelmarty.github.io/UniDecDocumentation/.

My goal is that this documentation will allow you to utilize the power of the UniDec python engine for scripting data analysis routines and performing custom analysis. Also, it will allow you to add new modules to the UniDec GUI.

I'm still working on documenting some of the windows and extensions (including MetaUniDec and C code), 
but the core features should be here.

##UniDec Architecture


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


##Getting Started


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

v. 2.2.1

Public Release Candidate

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








