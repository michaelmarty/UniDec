UniDec: Universal Deconvolution of Mass and Ion Mobility Spectra 
=================================================================

UniDec is a Bayesian deconvolution program for deconvolution of mass spectra and ion mobility-mass spectra.

It was orignally published in: [M. T. Marty, A. J. Baldwin, E. G. Marklund, G. K. A. Hochberg, J. L. P. Benesch, C. V. Robinson, Anal. Chem. 2015, 87, 4370-4376.]
(http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00140)

Detailed descriptions of the algorithm are provided in the paper. Please cite us if you use UniDec in your research.

UniDec may be downloaded from [unidec.chem.ox.ac.uk](http://unidec.chem.ox.ac.uk/).

Please contact michael.marty@chem.ox.ac.uk for questions, suggestions, or with any bugs.

## Downloading the Binaries

As described below, the Python code presented here relies on two critical binaries, UniDec.exe and UniDecIM.exe.

These binaries are available for download at [unidec.chem.ox.ac.uk](http://unidec.chem.ox.ac.uk/). 

These binaries should be deposited in the /unidec_bin directory once downloaded. 
 
If you want to convert Waters .Raw files, you will also need to add cdt.dll (for IM-MS) and MassLynxRaw.dll (for MS) to the same directory.
These files can be found [here](http://www.waters.com/waters/supportList.htm?cid=511442&locale=en_GB&filter=documenttype|DWNL&locale=en_GB) with support numbers DWNL134825112 and DWNL134815627.
Unfortunately, the libraries provided by Waters do not work on certain instruments, such as the G2-Si. 

##UniDec Documentation

This documentation is for the Python engine and GUI.

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







