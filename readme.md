UniDec: Universal Deconvolution of Mass and Ion Mobility Spectra 
=================================================================

UniDec is a Bayesian deconvolution program for deconvolution of mass spectra and ion mobility-mass spectra.

It was orignally published in: [M. T. Marty, A. J. Baldwin, E. G. Marklund, G. K. A. Hochberg, J. L. P. Benesch, C. V. Robinson, Anal. Chem. 2015, 87, 4370-4376.]
(http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00140)

Detailed descriptions of the algorithm are provided in the paper. Please cite us if you use UniDec in your research.

UniDec may be downloaded from [unidec.chem.ox.ac.uk](http://unidec.chem.ox.ac.uk/).

Please contact mtmarty@email.arizona.edu for questions, suggestions, or with any bugs.

## Installing

### Python

There are several python libraries that UniDec will depend on. 

matplotlib
numpy
scipy
wxpython
natsort
twython
pymzml
networkx

With the exception of wxpython, which can be installed from the [web](http://wxpython.org/), all of these can be installed from the command line with (for example):
    
    pip install natsort

Note: I would highly recommend setting up 64-bit Python as the default. MS data works fine with 32-bit, but IM-MS data is prone to crash the memory. If you are getting memory errors, the first thing to try is to upgrade the bit level to 64.

### Downloading the Binaries

As described below, the Python code presented here relies on two critical binaries, UniDec.exe and UniDecIM.exe.

These binaries are available for download at [unidec.chem.ox.ac.uk](http://unidec.chem.ox.ac.uk/). 

These binaries should be deposited in the /unidec_bin directory once downloaded. 
 
If you want to convert Waters .Raw files, you will also need to add cdt.dll (for IM-MS) and MassLynxRaw.dll (for MS) to the same directory. These files can be found [here](http://www.waters.com/waters/supportList.htm?cid=511442&locale=en_GB&filter=documenttype|DWNL&locale=en_GB) with support numbers DWNL134825112 and DWNL134815627. See below for more specifics. 

I have binaries built for Mac and Linux as well. They are a bit slower than the Windows version because they are compiled with gcc rather than the Intel C Compiler, but they are perfectly functional and still pretty darn fast. I can send these to you on request.

## Compatible File Types

UniDec is built to open .txt files using [numpy.loadtxt](http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.loadtxt.html). 

For MS data, it opens a two-column either a tab or space delimited list of m/z and intensity values.

For IM-MS, it will open a three-column tab or space delminated list of m/z, arrival time (or bin), and intensity values. Sparse matrices are fine for IM-MS. 

Recent versions are compatible with a text header at the beginning of the file. It will skip lines until it reaches the start of the data.

For Water's .raw files, UniDec is bundled with converters (CDCReader.exe and rawreader.exe) to convert the data to .txt. It will compress the retention time dimension into a single spectrum. A single file can be opened directly, or mulitiple files can be converted using Tools > Simple Batch Process Raw to Txt. For a fancier conversion such as extracting specific functions or scans, try Tools > Raw to Txt Conversion Wizard.

Water's converters will need [MassLynxRaw.dll](https://www.waters.com/waters/support.htm?locale=en_GB&lid=134815627&cid=511442&type=DWNL) (32-bit, for MS) and/or [cdt.dll](https://www.waters.com/waters/support.htm?locale=en_GB&lid=134825112&cid=511442&type=DWNL) (64-bit, for IM-MS) in the same directory as the converter executables (the unidec_bin folder or the top directory).

For Thermo .raw files, we do not currently have the conversion libraries built in (an enterprising developer should take this on...). 
The best way to load these files is to convert to mzML using [Proteowizard](http://proteowizard.sourceforge.net/). UniDec will open mzML file as if they are a text file.
We utilize [pymzML](http://pymzml.github.io/intro.html#general-information) for this. Cite them at: *Bald, T., Barth, J., Niehues, A., Specht, M., Hippler, M., and Fufezan, C. (2012) pymzML - Python module for high throughput bioinformatics on mass spectrometry data, Bioinformatics, doi: 10.1093/bioinformatics/bts066.*

##UniDec Documentation

Documentation is for the Python engine and GUI and can be found at http://michaelmarty.github.io/UniDecDocumentation/.

My goal is that this documentation will allow you to utilize the power of the UniDec python engine for scripting data analysis routines and performing custom analysis. Also, it will allow you to add new modules to the UniDec GUI.

I'm still working on documenting some of the windows and extensions, but the core features should be here.

The C code is proprietary and not open source. The copyright is owned by the University of Oxford. 
If you are interested in the C source code, please contact me. 
It should be possible to share with appropriate licensing or research agreements.

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

##Licensing

UniDec is free for noncommercial use under the [UniDec Academic License](http://unidec.chem.ox.ac.uk/12116_UniDec_Academic%20Use%20Licence.pdf)

Commercial licensing is available. Details are provided at the bottom of the Academic License.

The Python GUI and engine are licensed under the [GNU Public License (v.3)](http://www.gnu.org/licenses/gpl-3.0.en.html).

## Change Log

v. 1.0.11

Thanks to Tim for bug fixes on CDCReader.exe (the binary that converts Water's IMMS files to text files).

Updates to UniDecIM.exe to change how the twaveflag is determined and handled. Will now accept multiple possible twaveflags to allow for alternative calibration strategies in the future. Contact me with your function of interest.

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








