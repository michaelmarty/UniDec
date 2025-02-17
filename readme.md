UniDec: Universal Deconvolution of Mass and Ion Mobility Spectra 
=================================================================

UniDec's Mission: Making it easier to do more with your data.

UniDec is a Bayesian deconvolution program for deconvolution of mass spectra and ion mobility-mass spectra.

It was originally published in: [M. T. Marty, A. J. Baldwin, E. G. Marklund, G. K. A. Hochberg, J. L. P. Benesch, C. V. Robinson, Anal. Chem. 2015, 87, 4370-4376.](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00140)

Detailed descriptions of the algorithm are provided in the paper and subsequent papers. Please cite us if you use UniDec in your research.

Please contact mtmarty@email.arizona.edu for questions, suggestions, or with any bugs.

## Installation

UniDec may be downloaded from [https://github.com/michaelmarty/UniDec/releases](https://github.com/michaelmarty/UniDec/releases).

This compiled version is compatible with 64-bit Windows. It a portable binary, so it does not need a conventional installation.
Just unzip the folder, put it somewhere convenient, if you are still having issues 
run the INSTALLER.cmd as administrator,(only have to do this once)
then click the GUI_UniDec.exe file in the folder to launch.

To use the PDF report generator, install [MikTex](https://miktex.org) and select install packages automatically from the installation options. You may need to add this to your system path so that PDFLatex is found from the command line. 

## Tutorial

You can watch a video tutorial on how to use UniDec here: [https://www.youtube.com/watch?v=e33JxgY6CJY](https://www.youtube.com/watch?v=e33JxgY6CJY).

## Licensing

UniDec is distributed under a completely open source license. Our hope is that this allows UniDec to be
more widely used. If you are interested in including UniDec in another academic or commercial software distribution, 
you are welcome to email mtmarty@arizona.edu for more information. 

UniDec source code and compiled binaries are released under a modified BSD License as described below. Note, we ask
that you cite us in any publications. Quantitative citation metrics will help grant applications to support future development.

By downloading UniDec, you are agreeing to the UniDec and any third party licensing agreements and all terms therein. 

### UniDec License:

Copyright (c) 2016, University of Oxford
              2017-2023, University of Arizona
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

Intel libraries from the oneAPI are distributed under licenses found in the bin folder. 

RDKit libraries are gratefully used for LipiDec and are distributed under the CC-SA license. Please [cite them when using these features](https://www.rdkit.org/docs/Overview.html#citing-the-rdkit). 

FFTW libraries are distributed under the GPL license, which is included alongside the relevant files.

## UniDec Compatible File Types

UniDec is built to open .txt files using [numpy.loadtxt](http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.loadtxt.html). 

For MS data, it opens a two-column either a tab or space delimited list of m/z and intensity values.

For IM-MS, it will open a three-column tab or space delimited list of m/z, arrival time (or bin), and intensity values. Sparse matrices are fine for IM-MS. 

It is compatible with a text header at the beginning of the file. It will skip lines until it reaches the start of the data.

For Water's .raw files, UniDec is bundled with converters (CDCReader.exe) to 
convert the data to .txt. It will compress the retention time dimension into a single spectrum. 
A single file can be opened directly, or multiple files can be converted using 
Tools > Simple Batch Process Raw to Txt. For a fancier conversion such as extracting specific functions or scans, 
try Tools > Raw to Txt Conversion Wizard. Note: rawreader.exe has been replaced with the MassLynxSDK 4.5 Python library. Water's converters will need MassLynxRaw.dll and/or cdt.dll in the same directory as the converter executables (the unidec/bin folder or the top directory). You can find these at: [https://interface.waters.com/masslynx/developers-area/sdks/](https://interface.waters.com/masslynx/developers-area/sdks/) if they aren't already there. 

Agilent .d files can be read natively on Windows. 

Thermo .raw files should be able to be opened natively on Windows. Thermo DLLs are included bundled with UniDec.

Finally, many vendor formats can be converted mzML using [Proteowizard](http://proteowizard.sourceforge.net/). UniDec will open mzML file as if they are a text file, and this format should be cross platform.
We use [pymzML](http://pymzml.github.io/intro.html#general-information) for this. Please [cite them](https://www.ncbi.nlm.nih.gov/pubmed/22302572).

If you are a fan of mzXML, we recently added mzXML support courtesy of pyteomics (please cite: Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6 and Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717).

## MetaUniDec and UniChrom File Types and Importing

With MetaUniDec and UniChrom, everything is stored in a single HDF5 files. 
The HDF5 Import Wizard allows you to import a range of different file types directly into a single HDF5 file for MetaUniDec.
Thermo RAW and mzML files are supported fully, which means that either the scan or time range can be specified.
Text and Waters RAW files are supported for file import. Text files must be a single m/z spectrum.
Waters RAW files will have all scans summed into a single m/z spectrum upon import in MetaUniDec. 
The File>Waters Conversion Wizard tool allows specific scans to be converted into text files for importing.

In addition to the Import Wizard, there are several Manual File options, which will allow you to create blank HDF5 
(New File) and load data into it (Add File). Note, Add Data will sum all scans together, and Waters data is not supported.
You can select multiple files at once here. 
You can also just copy the data from XCalibur or MassLynx and then use Add Data From Clipboard. 

There are a few automated tools to parse single chromatograms directly into HDF5 files if you have all the data chromatograms 
with predictable scans or times. You can batch process multiple files at once. 

However, if you want to look at chromatography data directly, we now recommend UniChrom, which will plot the TIC and allow you to select specific scan ranges either manually or automatically.


## Installing the Source Code

Most users will likely just want to run the compiled version. For those advanced users who have experience with Python,
we have provided the source code for the GUI and API. For more information, check out this [walkthrough](https://github.com/michaelmarty/UniDec/wiki/Installing-the-UniDec-Source-Code). Specific package requirements are outlined in the setup.py file.

However, an experimental distribution is provide on PyPI. Using:

    pip install unidec

should install UniDec on your computer and install any required dependencies. Try this out first and see if it works.

### The UniDec Binaries

As described below, the Python code presented here relies on one critical binary, UniDec.exe. The binary should be in the /unidec_bin directory. 

If you are interested in building the binary or modifying the source code, the code and Visual Studio project files
are in unidec_src/UniDec. It is currently configured for Visual Studio Community 2022 with HDF5 1.12.2 and the Intel oneAPI compiler.
It can be easily compiled with other compilers but will show a significant performance loss without the Intel Compilers.
 
If you are interested in binaries for Mac and Linux, I would recommend building them yourself using the scripts in the unidec/src directory. UniDec compiles easily and works fine with these operating systems. I have provided some linux and mac version for you to test, but no promises they will work. Alternatively, you can use the [Docker container](https://hub.docker.com/r/michaeltmarty/unidec). 

## UniDec Documentation

Documentation is for the Python engine and GUI and can be found at http://michaelmarty.github.io/UniDecDocumentation/. Sorry, it's no great documentation, and I'm way behind on completing it. I'm still working on documenting some of the windows and extensions (including MetaUniDec and C code), but the core features should be here.

My goal is that this documentation will allow you to use the power of the UniDec python engine for scripting data analysis routines and performing custom analysis. Also, it will allow you to add new modules to the UniDec GUI.

## UniDec Architecture


UniDec is bilingual. The core of the algorithm is written in C and compiled as a binary.
It can be run independently as a command line program fed by a configuration file.

The Python engine and GUI serve as a very extensive wrapper for the C core. 

The engine (engine.py) can be operated by independently of the GUI. This allows scripting of UniDec analysis for more complex and high-throughput analysis than is possible with the GUI.
The engine contains three major subclasses, a config, data, and peaks object.

The GUI is organized with a Model-Presenter-View architecture.
The main App is the presenter (GUniDec.py).
The presenter contains a model (the UniDec engine at engine.py) and a view (mainwindow.py). 
The presenter coordinates transfer between the GUI and the engine.

MetaUniDec has a similar structure with a presenter (MetaUniDec.py), engine (mudeng.py), and view
 (mudview.py). However, unlike conventional UniDec, MetaUniDec includes a number of additional features to process 
 and analyze data. It relies less heavily on the Python API.

## Getting Started with UniDec in Python

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

You can also run files for simple deconvolution directly from the command line with:

    python -m unidec -f "C:\data\test.txt"

The main GUI class is Launcher. You can launch the Launcher from the command line with:

    python -m unidec.Launcher

Of course, using the pre-compiled version means you don't need to know Python at all and can just click the GUI_UniDec.exe icon to launch the program and get started. However you choose to do it, happy deconvolving!

## Change Log

v.8.0.1

Improved CD-MS defaults.

Added scan parsing to UniChrom. 

Added checkboxes to UPP to allow for more control over the HTML output. Moved some output controls down to better organize the window.

Fixed bug with dialog windows failing to return properly. 

Fixed a weird mzML bug.

v.8.0.0

Added **IsoDec** for isotopic deconvolution. This is a new window that allows you to deconvolve isotopic distributions. It is built using a neural net engine but has a similar interface to conventional UniDec.

New IsoGen neural net for generating isotopic distributions. Used by IsoDec but also available in Python code and via DLL.

Complete rewrite and restructuring of Importing code into UniDecImporter library. This should help harmonize and centralize the MS file reading functions. Removed some prior reliance on multiplierz for Agilent data.

To simplify the import functions, the UniDec GUI has changed slightly such that you need to explicity switch it to ion mobility mode (see Advanced menu) before opening ion mobility files.

Various other under the hood changes to accomodate new code, and some bug fixes.

Fixed label for DAR mode in UPP help files. 

v. 7.0.3

Added check box in UPP to allow all PNG images rather than SVG images, which can improve file sizes.

Added "Config Data Subbuff" to UPP to allow background subtraction with a curved background in the spreadsheet format. Added "Config Data Binevery", "Config Data Lintype", and "Config Data Smoothing" to provide access to other data processing parameters in the spreadsheet format.

v.7.0.2

Improvements to UPP outputs. Fixed issue with reports getting overwritten if the same file is used in multiple rows with UPP. 

Added new options in MetaUniDec for Average Charge State Plotting. Added exports of text files for fits and fitdata.

Added new pre-import intensity filter for UCD, which should help with memory crashes on large data files.

Fixed bug opening mzML files with new build.

Fixed bug with smashflag on UPP.

Fixed a few display and interface bugs on UCCD and UCD.

Fixed Waters Data Conversion Wizard issue with Trapping Collision Energy (eV). Stopped it from closing by default when the conversion is initiated.

Fixed bug with KD fitting on DataCollector. Fixed bug with UltraMeta. 

Fixed bugs with UCD Thermo imports and with Thermo header information for resolution.

v.7.0.1

Continued development and tweaks to UniChromCD. Bug fixes and GUI improvements. Added multiplex spectrum cropping. New HT sequences.

Switched default font to Arial after years of fighting Illustrator on Deja Vu Sans.

Bug fixes for file opening on Agilent and mzML.

v.7.0.0

Added new module for time domain CD-MS data analysis, UniChromCD. This includes some cool new plotting and analysis features. Most importantly, it also includes our new method for Hadamard Transform analysis of CD-MS spectra. You can use the tools for regular LC-CD-MS or go crazy with HT-LC-CD-MS. More info coming soon on this front...:)

UPP added NumPeaks output on all to help with catching bad files in the results. Also, fixed bugs to allow files with bad data to proceed relatively unscathed. It should be more error tolerant now for large runs. 

Added Ctrl+0 as shortcut to open the most recent file.

Added experimental feature for check box on main UniDec panel for whether to normalize the peak threshold. When on, it will behave as previously where the peak thershold is a ratio of the max intensity. When off, it will be a fixed intensity. Try switching off "Publication Mode" to view the plot with intensity axes. Note, this interacts with the Peak Normalization parameter directly above it. Max will normalize everything to 100. Total can get weird. None is the safest bet. 

Blocked mouse wheel events on drop down selections. This should fix a lot of accidental switches to parameters. 

Removed or made optional some dependencies to streamline installation. Removed special report type. Updates for Python 3.12 and general library updates.

Fixed deep bug with integration transforms. Improvements and refactoring to support UniChromCD. Fixed execution bug on servers. Fixed other bugs.

v.6.0.5

Merged into v.7

Minor fixes to update compatibility for Python 3.11 and 3.12. 

Fixed common crash bug when isotope mode is on. Added warnings for if isotope mode is on.

v.6.0.4

Added "Notes" option to UPP to add notes to the bottom of the report. Added peak and data normalization options to UPP spreadsheet inputs. See help file. 

Added "PyroGlu" option in UPP to automatically account for pyroglutamate formation on sequences that could have it. See help file for more information.

Printing HTML reports to PDF now looks better, with page breaks and background colors.

Command line now runs DScore with peak picking, so DScores will now appear in UPP.

Added del_columns custom field in batch.py for those that want to edit the report fields.

Fixed bug with CDCReader path in ImportWizard. Fixed weird bug with compressed mzML files. Fixed bug with save figure as and color bars.

v.6.0.3

Added global peak output to UPP. Updated documentation and help for UPP to show how to write your own workflows.

Added cosine similarity to Image Plotter.

Added "Config Smash File" input to UPP.

Added DScore to peakparam.dat output.

Fixed bug with Waters time range on file open.

Fixed bug with peaks shapes not registering. Note, this fix will break compatiblity for the charge peak shape
parameter on older config files for UCD. Just change zpsfun to zpsfn on the config file to rescue.

Fixed bug with PDF saving not working.

v.6.0.2

**Added a new Drug-to-Antibody Ratio (DAR) calculation mode for UPP.** See help file for the new keywords. 

Other improvements to UPP:
* Added color to the columns on UPP.
* Added ability to do peak integration.
* Added low and high charge config options.
* Added buttons to open results in Excel or other spreadsheet software. 
* Added global HTML report. 
* Improved file naming of outputs.
* Added progress bar on the bottom.

Fixed bug with HTML float sorting as text. Fixed major bugs with file imports.

v.6.0.1

**Added Noise Removal Tool** to UniDec and UniDecCD. This allows you to select specific noise regions (common in UHMR data) and set them to 0. Try Tools > Select Noise Peaks and the checking "Remove Noise Peaks" in the Advanced Data Processing Parameters.

UPP Changes: 
* Added "Config m/z Peak FWHM", "Config m/z Peak Shape", and "Config Sample Mass Every" parameters. If "Config m/z Peak FWHM" is specified as a float or integer, it will use that and skip the automatic peak width determination.
* Enabled relative paths. 
* Added "All" or "None" keywords on "Apply Fixed Mods" and "Disulfides Oxidized." 
* Added "Config File" keyword to load a specific config file.
* More tolerant of nans and different file extensions
* Added ability to select files to load to UPP.
* Added "Global Fixed Mod" to apply the same float value to each pair.
* Special Bispecific Antibody Pairing Calculations (see Help File)
* Fixed issues bugs with use_converted, clear all, and resetting the config.

Fixed major bugs with 2D zooming and plotting. 

v.6.0.0

Major changes to code structure to enable PyPI distribution. Changed to more conventional Python package structure with unidec top directory. Renamed unidec_modules to modules, unidec_src to src, and unidec_bin to bin. Renamed some top level windows and moved all to main folder. For example, unidectools.py has been moved from modules to the top level and renamed as tools.py, and unidec.py has been renamed as engine.py. The good news is that you can now do `pip install unidec` and `python -m unidec.Launcher` to run UniDec from the command line.


Improved command line tools. Added -o option for file output in unidec. Included new entry points and smoothed out command line interface. Now you can use unidec to lanuch python -m unidec in the scripts folder. Similarly, gunidec will launch python -m unidec.Launcher. If you have Python and the scripts folder configured in your system path, command line tools will be much easier to use. 

Added a new mode to the Oligomer and Mass Tools, Site Mode. Here, you can specify specific binding sites and which potential species could be bound in that site. Any number of rows (species) or columns (sites) can be added as long as your memory can handle the possible combinations. Use 0s and 1s to incidate whether each species can be bound at each site. You can copy/paste from excel or open a CSV or XLSX file. 

Added the first build of the **UniDec Processing Pipeline** tool to help with batch processing and analysis for large sample sets. More details to come.

Added new experimental HTML report generator. Found in File > Save Figure Presets. Test this out and let me know what you'd like to see in it.

Added a color plot option to the right click of the peak panel. Check it out!

Added FWHM and error to the Copy All Full right click option in the peak panel.

**Changed default behavior for zooming such that it does not zoom out when plots are replotted.** This is a major change that I think will be more intuitive. If you want to zoom out, you can still click once on the plot as usual. 

Removed the ability to manually integrate peaks. Now, only auto integration is allowed. If you want to have this added back in, let me know, but it seemed to be causing more trouble than it was worth.

Various bug fixes including fixing EPS export and bugs on Export2D. Fixed import of tab delimited files. Resized manual file menu. Improved window stretching in oligomer and mass tools. Broken plot peaks on UniDecCD. Bug on selecting alternative matches not propagating. Updated background color to blue on alternative matches to avoid confusion with RYG scheme. Fixed labeling bug with isolated match names. 

v.5.2.1

UniDec and UniChrom now try to read the polarity from the data and set the adduct appropriately. Check that this is working!

Added right click for zoom on UniChrom TIC.

Added ability to label names from peak panel.

Added drag and drop for batch processing in UniChrom. 

Significantly faster reading of some mzML files. 

Fixes to reading in Agilent Files. May no longer need adminstrator access but not sure. Also, may work with Bruker TIMS data, but also not sure about this.

Fixes with reading some config files. Fixes with npz file reading.

v.5.2.0

Major speed up to Plot Peaks command that simulates the isolated charge state distributions! This function is now built into the same C code used by UniDec, and is called by the UniDec.exe <file>_conf.dat -conv. It produces a binary file output that is read in and plotted. 

Improvements to matching speed and API. Added right click option to send item from oligomer list to common matches list. Added a bunch of glycans and stuff to the common masses table.

Improvements to MetaUniDec. Create a merged template rather than just using the first spectrum. Improvements to speed by keeping file open for the whole session. Reduction in print outputs.

Improvements to Imaging viewer, including bug fixes.

Added ability to calculate consecutive differences in peak masses.

Added continuous plot button for UCD m/z vs. mass. Also, fixed bug to improve plotting.

Improved copy/paste from peak list.

Added ability to read text file exported from Chromeleon that have a bunch of commas in the numbers.

Added mzXML support courtesy of pyteomics (please cite: Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6 and Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717)

Added Docker container and updated for command line usage of unidec.py.

Fixed bug with the Peak Width Tools in UCD. Fixed labeling issue on UniChrom. Added UCD as command launch item.

v.5.1.1

Added Smashing (Ctrl+Double Right Click) to UCD.

Fixed deep bugs with HDF5 in C code by going back to 1.12.0. Reverted C code back to without _s functions for Mac and Linux compatibility. Fixed bugs with report generator.

v.5.1.0

Added Intensity Threshold to MetaUniDec. Updated code throughout to allow more tolerance for empty data sets below threshold. Streamlined some print commands.

Ported UniDecCD from Python to C, which speeds things up. The point smoothing is slightly different, but it should otherwise behave very similarly.

Added ability to label the peak areas on the deconvolved mass plot. Check out the right click menu in the peak list.

Fixed a bunch of warnings when compiling on VS2022. Updated to _s libaries for a lot of stdlib stuff. If this breaks compiling on Mac and Linux, let me know.

Added experimental MS Imaging tools.

Added tool for identifying alternative matches on the Oligomer and Mass Tools and for selecting alternate matches from the list. Explore right clicks on the match panel to test this out.

Speeding up a few things in Meta.

Fixed bug with multiple file dialog. Fixed bug with biopolymer calculator. Fixed bug with peak list. Fixed bug where it was crashing for people in certain regions (update to wxpython).

v.5.0.5

**Upgraded build** using the Intel oneAPI and Visual Studio 2022. Giving about a 2x speed boost!

Added new Gaussian fitting experimental features to DataCollector.

Fixed bug in build with PDF export.

Fixed bug with manual assignments.

Fixed CDMS example data showing up in other windows.

v.5.0.4

Added Mass Limits to Mass Defect Tools.

Added Normalize Data option in UniDecCD.

Allowed Waterfall plots to show text data on the y-axis.

Fixed bug in Thermo Importer. Fixed bugs with a few parameters in UCD.

Added ability to open .dmt and .i2ms files from STORI analysis.

Added experimental MassQL feature to select peaks that match certain queries. 

Adding experimental High Throughput Screen features in MetaUniDec.

Added ability to open npz files, which are much faster than txt and should help in custom scripts.

v.5.0.3

Added m/z vs. mass plot for UniDecCD using a button in the Additional Plotting Parameters. Thanks to Sean Cleary for inspiration here.

A few new experimental features in my attempt at real time deconvolution:

* Added refresh on UCD.

* Added autorefresh tool on UniDec and UniChrom. If you open a file that is currently collecting data, you can do real time deconvolution as the data is collected.

* A command line argument will automatically launch to UniChrom and open the file argument that is provided. This means you can now set UniDec as the default app in Windows to open Raw files and other data format. If you do this, when you click the icon on the Thermo Exactive Tune Software, it will open UniDec as the default app. Still a work in progress, but kind of fun to play with. 

Added new features to Extract 2D window. 

Finally figured out how to fix the bug for wxPython 4.1, but it involves modifying the source code. Contact me if you are interested in running it yourself. 4.1 will be reflected in the build.

v.5.0.2

Added Batch Processing (via Tools menu or drag and drop) to UniDecCD.

Fixed bugs with STORI calibration in UCD.

Added export of text file for composite spectrum.

v.5.0.1

**Added Beta Support for STORI Folders** in UCD. Select the directory of CSV files under Tools. It will concatenate these and convert them into an npz file that you can open directly.

**Added Beta Support for Agilent Drift Tube IM-MS**. To use this, convert your data to mzML using MSConvert and select the "Combine Ion Mobility Scans" option. We recommend using the package as gzip as well and setting the extension to ".mzML.gz". Next, open UniDec and switch to Ion Mobility Mode under the Advanced Tab. Open the mzML file (drag and drop will work). Adjust the parameters in data processing to get it to look nice. The pusher should be set to 0. Set the voltage, temperature, and pressure to 0 in the Ion Mobility Parameters (with Linear Cell selected). Enger the tfix as the Dead Time and the Beta parameter below in the Drift Cell Length/Beta box. Adjust your parameters and hit deconvolve.

**Added Support for SLIM TWIMS IM-MS**. Select SLIM poly3 or poly2 under the T-Wave Calibration Type. Calibration parameter 1 is the constant term with parameter 2 as the linear term and so on. 

Added ability to export m/z values for peaks as _mzpeakdata.dat. Added headers to peakparams export text file.

Fixed bugs with UltraMeta.

v.5.0.0

**Added UCD: UniDec for Charge Detection-MS**. This major new window extends UniDec to CD-MS data. It builds on the existing UniDec GUI but uses a new deconvolution engine. The engine is written in Python and has GPU-acceleration available for anyone with CUDA 11.2 installed. 

**Added mzML ion mobility support**. Using MSConvert, select "Combine ion mobility scans". Then, launch UniDec and switch to IM mode by clicking Advanced > Switch to Ion Mobility Mode. UniDecIM should then be able to open ion mobility mzML files by drag and drop or File > Open File. It will create a text file next to the mzML file that you can open in the future. We recommend using the "Compress when converting to .txt" option to help speed up data processing by binning the data at this stage using the "Bin Every" parameter. 

**Improved DoubleDec** code that is integrated into the C code and more tolerant of different types of input.

Began adding experimental SLIM IM-MS calibration functions. Not fully implemented.

Added ability to open CSV files exported from Thermo Freestyle.

Cleaned up and refactored some code to fit UCD. Bug fixes.

v.4.4.1

Added ability to open UniDec GUI from UniChrom selection and ability to right click a spectrum and add it to selection on UniChrom. 

Added time limit option to UniChrom.

Added Subtract and Divide window to help visualize average monomer masses.

Added commandline options to the Python code to allow -f "file.txt" to automatically launch the file with the program. Using -u/--unidec, -m/--meta, or -c/--chrom will launch UniDec, MetaUniDec, and UniChrom respectively with the launcher.

Implemented DoubleDec in C code (Jack Liu) rather than Python.

Edited Waters Raw to Text Conversion Wizard to allow specification of extraction times rather than just scans.

Fixed bug with FPOP. Fixed bug with charge extraction in C code. Fixed bug with single scan Thermo files. Fixed bug with remove duplicates. 

Added a command line warning when data might be approaching the memory limit.

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






