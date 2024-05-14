import wx
import wx.html
import os


class HelpDlg(wx.Frame):
    def __init__(self, num, *args, **kw):
        super().__init__(*args, **kw)
        pathtofile = os.path.dirname(os.path.abspath(__file__))
        self.imagepath = os.path.join(pathtofile, "images")
        # print(pathtofile)
        # print(self.imagepath)
        if num == 1:
            self.get_started()
        elif num == 2:
            self.import_window()
        elif num == 3:
            self.plot_window()
        elif num == 4:
            self.peak_window()
        elif num == 5:
            self.data_processing()
        elif num == 6:
            self.unidec_parameters()
        elif num == 7:
            self.additional_filters()
        elif num == 8:
            self.peak_selection()
        elif num == 9:
            self.additional_plotting()
        elif num == 10:
            self.auto_import_help()
        elif num == 11:
            self.presets_help()
        elif num == 12:
            self.batch_help()
        elif num == 13:
            self.peak_width_tool_help()
        elif num == 14:
            self.oligomer_help()
        elif num == 15:
            self.auto_match_help()
        elif num == 16:
            self.animate_help()
        elif num == 17:
            self.autocorr_help()
        elif num == 18:
            self.fft_help()
        elif num == 19:
            self.baseline_help()
        else:
            html = wx.html.HtmlWindow(self)
            html.SetPage("<html><body>You shouldn't see this!!! ERROR!!!!</body></html>")

    def get_started(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(600, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Welcome to MetaUniDec!</h1>"
            "<p>UniDec is a fast, robust, and flexible software package "
            "used for the deconvolution of mass spectra and ion mobility-mass spectra.</p>"
            "<h2>MetaUniDec File Types and Importing</h2>"
            "<p>With MetaUniDec, everything is stored in a single HDF5 files."
            "The HDF5 Import Wizard allows you to import a range of different file types directly into a single HDF5 file."
            "Thermo RAW and mzML files are supported fully, which means that either the scan or time range can be specified."
            "Text and Waters RAW files are supported for file import. Text files must be a single m/z spectrum."
            "Waters RAW files will have all scans summed into a single m/z spectrum upon import."
            "The File>Waters Conversion Wizard tool allows specific scans to be converted into text files for importing.<p>"

            "<p>In addition to the Import Wizard, there are several Manual File options, which will allow you to create blank HDF5"
            "(New File) and load data into it (Add File). Note, Add Data will sum all scans together, and Waters data is not supported."
            "You can select multiple files at once here."
            "You can also just copy the data from XCalibur or MassLynx and then use Add Data From Clipboard.<p>"

            "<p>There are a few automated tools to parse single chromatograms directly into HDF5 files if you have all the data chromatograms"
            "with predictable scans or times. You can batch process multiple files at once."
            "Only mzML and Thermo RAW files are supported for automated chromatogram import.<p>"

            "<p>Note: For import Thermo RAW Data, you need msfilereader (and multiplierz if not running the .exe).</p>"

            "<h2>Loading Data</h2>"
            "<p>To open files with MetaUniDec, just drag and drop the HDF5 file into MetaUniDec, "
            "or go to File->Open File. You can also drag and drop and HDF5 file into the main window.</p>"
            "<h2>Analyzing data</h2>"
            "<table style=\"width:100%\"><tr><td>For basic data processing, deconvolution, and peak selection, simply "
            "hit the \"All\" button in the top right corner. For more rigorous analysis, see the other help guides.</td>"
            "<td><img src=\"" + self.imagepath + "/allButton.png\" alt=\"PNG Icon\"></td></table>"
                                                 "<h2>Saving Figures</h2>"
                                                 "<p>There are 3 ways to save your figures. To save an individual plot, middle-click the plot and a dialog "
                                                 "window will pop up that allows you to save it. The other methods save all available plots. "
                                                 " File->Save Figure As will allow you to select the directory "
                                                 "and file name header, along with the extension and image dimensions. File->Save Figure Presets will save "
                                                 "your figures as .pdf/.png/.pdf thumbnails to the default location, which is the location of the original "
                                                 ".HDF5 file. Ex: C:\\HDF5-location\\UniDec_Figures_and_Files\\Fig1.pdf<\p>"
                                                 "<h2>Importing Configs</h2>"
                                                 "<p>File->Load External Config will load the parameters from a previous run."
                                                 " This can either be a _conf.dat file from UniDec of an HDF5 file from MetaUniDec.  "
                                                 "Advanced->Reset Factory Default will restore all settings to the program defaults.<p>"
                                                 "</body></html>"
        )

    # TODO: Figure out what Make top does
    def import_window(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(500, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>The Spectra Table</h1><p>This page will teach you how to work with spectra within MetaUniDec, "
            "in the data table on the left side.</p>"
            "<h2>Right Click Options for the Spectra Table</h2>"
            "<img src=\"" + self.imagepath + "/rightClick.png\" alt=\"PNG Icon\">"
                                             "<table style=\"width:100%\">"
                                             "<tr><td>Ignore</td><td>Hides the sample on the graphs and table.</td></tr>"
                                             "<tr><td>Isolate</td><td>Singles out the selected sample on the graphs and table. Hides others.</td></tr>"
                                             "<tr><td>Repopulate</td><td>Brings back samples hidden by ignore/isolate to the graphs and table.</td></tr>"
                                             "<tr><td>Ignore</td><td>Hides the sample on the graphs and table.</td></tr>"
                                             "<tr><td>Analysis Tools</td><td>Does autocorrelation or FFT analysis on a single sample. "
                                             "See Help->Menu Bar->Analysis for a description of these tools.</td></tr>"
                                             "<tr><td>Change Color</td><td>Changes the color of a sample</td></tr>"
                                             "<tr><td>Make Top</td><td>Makes the spectrum the primary spectrum for peak fitting.</td></tr>"
                                             "<tr><td>Fill Down Variable 2</td><td>Changes all Variable 2 values to the selected sample's "
                                             "Variable 2 value.</td></tr>"
                                             "<tr><td>Delete</td><td>Deletes the sample. WARNING: This deletes the sample from the HDF5 file.</td></tr>"
                                             "<h2>Meta Data (Variables 1 and 2)</h2>"
                                             "<p>Variable 1 and 2 can be manually adjusted from the table. "
                                             "File->Variables allows the variables to be renamed, imported, and exported.<p>  "

                                             "</body></html>"
        )

    def plot_window(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Plot Window</h1><p>This page will teach you what each plot shows and how to manipulate the plots. "
            "Note: The plot window relies heavily on the controls described in Help->UniDec Controls.</p>"
            "<h2>Interacting with the plots</h2>"
            "<h4>Zooming in</h4><p>There are two ways to zoom in on the plots. First, by left clicking and "
            "dragging in a straight line, you will set the left and right bounds of the plot to the ends of the line. "
            "Second, by left clicking and dragging to create a rectangle, you will set the plot to the 4 edges of the "
            "rectangle. To undo the zoom, simply left click on the plot once.</p>"
            "<h4>Clicking on the plot (MS Data plot only!)</h4><p>On the unzoomed plot, single left click at two "
            "different spots on the plot and numbers will appear at those spots of the same color (indicate the "
            "charge), along with a number in the top right corner (indicates the mass). "
            "Holding the conrol key will allow the same actions on a zoomed window.</p><br><br>"
            "</p> Right clicking on this plot when zoomed will set the m/z range for data processing.</p>"
            "<h2>Plot Types</h2>"
            "<table style=\"width:100%\">"
            "<tr><td>MS Data</td><td>Plot of m/z vs. intensity. "
            "The highest peak shown will be set at intensity=1 if Normalize Data is checked and "
            "the rest of the peak's height will be normalized according to that intensity value.</td></tr>"
            "<tr><td>Charge Distributions</td><td>Plot of summed charge vs. intensity. This plot only shows data after running"
            " \"Run UniDec\".</td></tr>"
            "<tr><td>Mass Distributions</td><td>Plot of summed mass vs. intensity. This plot only shows data after running"
            " \"Run UniDec\". After running \"Peak Detection/Extraction\", a black line will appear representing "
            "average mass distribution. Peaks will be picked from this distribution and marked with the shape labels.</td></tr>"
            "<tr><td>Extracts Line Plot</td><td>Plot of variable 1 vs. intensity for the extracted peaks. This "
            "plot only shows data after running \"Peak Detection/Extraction\".</td></tr>"
            "<tr><td>Extracts Grid Plot</td><td>Heatmap of variable 1 vs. intensity (color represents intensity)"
            " for the extracted peaks. This plot only shows data after running \"Peak Detection/Extraction\".</td></tr>"
            "<tr><td>Bar Chart</td><td>Plot of intensity of each species (the peaks) at each variable 1. "
            "This plot only shows data after running \"Peak Detection/Extraction\".</td></tr>"
            "<tr><td>m/z Grid</td><td>Plot of m/z vs. variable 1 (color represents intensity). "
            "This plot only shows data after running \"Plot 2D Grids\".</td></tr>"
            "<tr><td>Mass vs. Charge</td><td>Plot of mass vs. variable 1(color represents intensity). "
            "This plot only shows data after running \"Plot 2D Grids\".</td></tr></table>"
            "</body></html>"
        )

    def peak_window(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Peak Window</h1><p>This page will teach you how to interact with the peak window.</p>"
            "<h2>Getting the window to appear</h2><p>The peak window is located to the right of the plot window. "
            "When you first load data into MetaUniDec, nothing will appear there. For data to appear in the peak "
            "window, you have to run \"All\" or \"Peak Detection/Extraction\".</p>"
            "<h2>Right Click Options</h2>"
            "<img src=\"" + self.imagepath + "/peakRightClick.png\" alt=\"PNG Icon\">"
                                             "<table style=\"width:100%\">"
                                             "<tr><td>Ignore</td><td>Hides the peak on the plots and table.</td></tr>"
                                             "<tr><td>Isolate</td><td>Singles out the selected peak on the plots and table. Hides others.</td></tr>"
                                             "<tr><td>Repopulate</td><td>Brings back peaks hidden by Ignore/Isolate.</td></tr>"
                                             "<tr><td>Label Charge States</td><td>On the MS Data plot, labels the different charge states for the "
                                             "selected sample with dashed-black lines.</td></tr>"
                                             "<tr><td>Display Differences</td><td>For the selected peak, gets the difference between it's mass and "
                                             "all of the other peak's masses. Will display these masses on Mass Distribution plot.</td></td>"
                                             "<tr><td>Color Select</td><td>Allows you to choose the color of the selected peak</td></tr></table>"
                                             "</body></html>"
        )

    def data_processing(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Data Processing</h1><p>Note: Hovering over the buttons/text fields will show helpful tooltips.</p>"
            "<table style=\"width:100%\">"
            "<tr><td>m/z Range</td><td>Limits the m/z range shown/analyzed.</td></tr>"
            "<tr><td>Baseline Subtraction</td><td>Uses the curved baseline option in UniDec. If the baseline intensity isn't 0, you can use this to "
            "subtract down to the desired baseline. Smaller numbers mean a more aggressive baseline subtraction. </td></tr>"
            "<tr><td>Bin Every</td><td>Averages <i>n</i> data points together. Can help smooth noisy data and reduces "
            "the number of data points. Since the speed of the UniDec algorithm depends mostly on # of data points and "
            "allowed charge states, binning can significantly increase the speed of the algorithm.</td></tr>"
            "<tr><td>Normalize Data</td><td>Selects whether your results will be normalized or not. Checking the box "
            "normalizes your data to a maximum intensity of one.</td></tr>"
            "<tr><td>Process Data Button</td><td>Processes the data according to how you filled out the above fields "
            "and displays the new plots when done.</td></tr></table>"
            "</body></html>"
        )

    def unidec_parameters(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>UniDec Parameters</h1><p>Note: Hovering over the buttons/text fields will show helpful tooltips.</p>"
            "<table style=\"width:100%\">"
            "<tr><td>Charge Range</td><td>Set the min and max charges for the deconvolution algorithm.</td></tr>"
            "<tr><td>Mass Range</td><td>Set the min and max masses for the deconvolution algorithm</td></tr>"
            "<tr><td>Sample Mass Every</td><td>Sets the distance between data points in the Mass Distribution plot."
            "<tr><td>Peak FWHM (Full Width Half Max)</td><td>FWHM of an expected peak. "
            "The Automatic Peak Width (Tools or Ctrl+W) will set this automatically. "
            "Tools->Peak Width Tools provides a fitting of the top spectrum.</td></tr>"
            "<tr><td>Peak Shape Function</td><td>Your expected peak shape. Split G/L is a Gaussian on the low m/z side "
            "and Lorentzian on the high m/z side to give a long tail at higher m/z.</td></tr>"
            "<tr><td>Run UniDec</td><td>Runs the deconvolution algorithm according to how you filled out the above "
            "fields. Note: after clicking this button, the Charge and Mass Distribution plots will show up.</td></tr>"
            "</table></body></html>"
        )

    def additional_filters(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Additional Filters/Restraints</h1><p>Advanced parameters for UniDec. Note: Hovering over the "
            "buttons/text fields will show helpful tooltips.</p>"
            "<table style=\"width:100%\">"
            "<tr><td>Charge Smooth Width</td><td>See tooltip. Almost always set to 1 unless using Mass Smooth or "
            "Isotope Mode instead, in which case it may be turned off if desired using 0.</td></tr>"
            "<tr><td>Mass Difference</td><td>Part of the mass smooth filter. Mass difference will incorporate "
            "neighboring species of known mass, such as neighboring oligomers, to help determine the charge.</td></tr>"
            "<tr><td>Mass Smooth Width</td><td>Best to use a width of 1 for on and 0 for off.</td></tr>"
            "<tr><td>Native Charge Offset Range</td><td>Limits charges to a windowed offset from the predicted "
            "native charge. Useful for eliminating extremely high or low charge states in complex samples. </td></tr>"
            "<tr><td>Isotope Mode</td><td>Uses Averagine model to project isotopic distributions, determine the charge "
            "state, and return the monoisotopic masses.</td></tr>"
            "<tr><td>Manual Mode</td><td>Forces m/z values within a defined window to be a defined charge. Opens "
            "Tools->Manual Assignment window.</td></tr>"
            "<tr><td>Charge Scaling Mode</td><td>See tooltip.</td></tr>"
            "<tr><td>Mass List Window</td><td>See tooltips. Opens Tools->Oligomer and Mass Tools window</td></tr>"
            "<tr><td>m/z to Mass Transformation</td><td>See tooltips.</td></tr>"
            "<tr><td>Maximum # of Iterations</td><td>See tooltip.</td></tr>"
            "<tr><td>Adduct Mass</td><td>See tooltip.</td></tr>"
            "</table></body></html>"
        )

    def peak_selection(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Peak Selection, Extraction, and Plotting</h1><p>Note: Hovering over the buttons/text fields will show "
            "helpful tooltips.</p>"
            "<table style=\"width:100%\">"
            "<tr><td>Picking Range</td><td>A peak will not be considered a peak unless it a maximum within"
            "+/- the picking range.</td></tr>"
            "<tr><td>Picking Threshold.</td><td>See tooltip.</td></tr>"
            "<tr><td>Peak Normalization</td><td>See tooltip.</td></tr>"
            "<tr><td>How to Extract Peaks</td><td>Select your method of extracting peaks. See tooltip for description "
            "of each type of extraction method.</td></tr>"
            "<tr><td>Extraction Window</td><td>See tooltip.</td></tr>"
            "<tr><td>Extract Normalization</td><td>Set how to normalize the extraction.</td></tr>"
            "<tr><td>Peak Detection/Extraction</td><td>Clicking this button will find and extract peaks according to "
            "the above values. THis button will make data appear in the Extracts Lines Plot, Extracts Grid Plot, and "
            "Bar Chart.</td></tr>"
            "<tr><td>Plot 2D Grids</td><td>Clicking this button will make the m/z Grid and Mass vs. Charge plots "
            "appear. Note: the \"All\" button by default does not execute this button, as this button can take a "
            "long time to run.</td></tr>"
            "</table></body></html>"
        )

    def additional_plotting(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Additional Plotting Parameters</h1><p>Note: Hovering over the buttons/text fields will show "
            "helpful tooltips.</p>"
            "<table style=\"width:100%\">"
            "<tr><td>2D Color Map</td><td>Selects color function for the m/z Grid and Mass vs. Charge plots.</td></tr>"
            "<tr><td>Peaks Color Map</td><td>Selects color function for the Extracts Line and Bar Chart plots, and "
            "changes the color of the peaks in the peak window.</td></tr>"
            "<tr><td>Discrete Plot</td><td>Uses discrete rectangles for the 2D plots."
            " This speeds up 2D plots but makes them less visually appealing.</td></tr>"
            "<tr><td>Publication Mode</td><td>Cleans up the plots and axes for publication.</td></tr>"
            "<tr><td>Reconvolved/Raw Results</td><td>See tooltip.</td></tr>"
            "<tr><td>Marker Threshold</td><td>Adds markers to peaks that are above this threshold.</td></tr>"
            "<tr><td>Species Separation</td><td>Set the amount of space separating different species in the MS Data, "
            "Charge Distributions, and Mass Distribution plots.</td></tr>"
            "<tr><td>Integration Range</td><td>Will integrate from peak-min to peak+max.</td></tr>"
            "<tr><td>Replot</td><td>Replot some of the plots. </td></tr>"
            "<tr><td>Plot Composite</td><td>Plots only the average mass distribution."
            "</td></tr></table></body></html>"
        )

    # TODO: Finish these

    def auto_import_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Auto Import Chromatograms</h1><p>This function will turn 1 or more chromatograms into HDF5 files. "
            "Chromatograms needed to be in .raw or .mzml format. Files need to be "
            "in the format filename_Ramp_startvoltage_endvoltage_voltagestep. Will create a HDF5 file with the same "
            "name.</p>"
            "<h3>Auto Import Chromatogram By Time</h3>"
            "<p>This function imports a chromatogram and separates it into datasets based on an inputted time step. "
            "Select a file. Select a time step to parse your file. Generally a voltage step of 20 corresponds to "
            "1.0 minutes, voltage step of 5 corresponds to 0.25 minutes, etc. Will create a HDF5 file with the same "
            "name.</p>"
            "<h3>Auto Import Chromatogram By Scans</h3>"
            "<p>This function compressed your chromatogram by averaging <i>x</i> number of scans into 1 data "
            "point, where <i>x</i> is the number you input. Select a file and enter a number to bin your file with.</p>"
            "<h3>Auto Import Multiple Chromatograms By Range of Times</h3>"
            "<p>This function takes multiple chromatograms and takes out an inputted time range from each file and "
            "combines them into 1 HDF5 file. Note: All files must have same time step/voltage step. Select files, "
            "enter the time step, followed by the start and end times in terms of the time step, and name the new HDF5 "
            "file.</p>"
            "<h3>Auto Import Multiple Chromatograms By Range of Scans</h3>"
            "Basically the same as Auto Import Multiple Chromatograms By Range of Times except you enter the start and "
            "end scans you want to import. Compiles all scans from 1 file into 1 dataset in the HDF5 file.</p>"
            "</body></html>"
        )

    def presets_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body><h1>Presets</h1>"
            "<p>The presets will change the Data Processing, UniDec Parameters, and Additional Filters/Restraints "
            "values to presets based on the type of mass spec being done. Additional adjustmets are likely necessary "
            "based on individual samples.</p>"
            "</body></html>"
        )

    def batch_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>There are several options for batch processing multiple HDF5 files with "
            "MetaUniDec. Batch Assign Configs simply applies the current config as defined by"
            " the controls on the right to all selected files. Batch Run will perform data processing "
            "and deconvolution on each file. Batch Extract will extract the peak information from each of the "
            "files provided the deconvolution has already been performed. "
            "Batch Assign/Run/Extract will do all of these at once.</body></html>"
        )

    def peak_width_tool_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>The Peak Width Tool allows fitting of the \"top\" spectrum, which by default is"
            " the first on the list of spectra. You can right click a spectrum and select \"Make Top\""
            " to make it the spectrum applied. Hitting Ok will update the main controls with the"
            " fit peak shape. Automatic Peak Width (Ctrl+W) will do this automatically.</body></html>"
        )

    def oligomer_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body><h1>Mass and Oligomer Tools</h1>"
            "<p>This tool allows you to match peaks to common species. Note: need to run Peak Detection/Extraction "
            "first. First, select the species you believe are in your "
            "mass spec data. If your species isn't in the table, you can import your own table with the Import from "
            "File button. File should be in a .csv file with the format Name, Mass, Type (column headers not "
            "required). Select your species, right click, and select Add to Oligomer Builder. In the Oligomer List, you"
            " can select the min/max # of oligomers for each species. You can also add your own oligomers with the "
            "Import from File button or Add Oligomer Species, which allows you to manually add the "
            "information about the oligomer. You can then click the Match to Isolated Oligomers button, which will try "
            "to match a peak to individual oligomers, or the Match to Mixed Oligomers, which will mix and match "
            "different oligomers together while trying to find a complementary peak. All peaks will have have the "
            "closest oligomer(s) compared to it, along with the error. Simulate Masses is unsupported in MetaUniDec.</p>"
            "</body></html>"
        )

    def auto_match_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Auto Match Peak</h1>"
            "<p>This function automatically matches the peaks to oligomers from the Oligomer and Mass tools. Will fill "
            "the oligomer list and try to match peaks to oligomer combinations.</p>"
            "</body></html>"
        )

    def animate_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>"
            "<h1>Animate</h1>"
            "<p>This function will show how the different plots change over Collision Voltage, timemid, etc. Can adjust"
            " how quickly plots change with the slider. Back and next buttons move to the previous/subsequent dataset "
            "respectively. Autoscale centers the plot.</p>"
            "</body></html>"
        )

    def autocorr_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>Performs an autocorrelation of the zero-charge mass spectrum to idenfity repeating mass units.</body></html>"
        )

    def fft_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>Performs two sequential Fast Fourier Transforms on the m/z spectrum to identify repeating units.</body></html>"
        )

    def baseline_help(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(400, 400))
        html = wx.html.HtmlWindow(self)
        html.SetPage(
            "<html><body>Under the Advanced menu, an experimental feature can be activated to determine "
            "the baseline automatically during the course of deconvolution."
            "</body></html>"
        )

        # TODO: add menu bar help stuff
