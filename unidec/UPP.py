import platform
import wx
import pandas as pd
from copy import deepcopy
from unidec.modules.isolated_packages.spreadsheet import *
from unidec.batch import UniDecBatchProcessor as BPEngine
from unidec.batch import *
from unidec.modules.html_writer import *
from unidec.GUniDec import UniDecApp


class HelpDlg(wx.Frame):
    def __init__(self, num=1, *args, **kw):
        super().__init__(*args, **kw)
        pathtofile = os.path.dirname(os.path.abspath(__file__))
        self.imagepath = os.path.join(pathtofile, "images")
        self.htmlstr = ""
        # print(pathtofile)
        # print(self.imagepath)
        if num == 1:
            self.help_frame()
        else:
            html = wx.html.HtmlWindow(self)
            html.SetPage("<html><body>You shouldn't see this!!! ERROR!!!!</body></html>")

    def help_frame(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title="Help", size=(600, 600))

        # Add a menu button to open in the web browser
        menu = wx.Menu()
        menu.Append(wx.ID_ANY, "Open in Browser")
        menu.Bind(wx.EVT_MENU, self.open_in_browser)
        menu_bar = wx.MenuBar()
        menu_bar.Append(menu, "Help")
        self.SetMenuBar(menu_bar)

        html = wx.html.HtmlWindow(self)

        html_str = gen_style_str()

        html_str += "<html><body>" \
                    "<header><h1>Overview</h1></header><p>" \
                    "Welcome to the UniDec Processing Pipeline (UPP)! " \
                    "This module is designed to help you process, deconvolve, " \
                    "and extract specific information from your data. " \
                    "Expanding on the batch processing features present in UniDec from the beginning, " \
                    "it is designed to interface with Excel/CSV files so that you can connect it with your workflows. " \
                    "</p>" \
                    "<h2>Basic Features</h2> <h3>Opening A File</h3><p>" \
                    "Although you can type everything directly into UPP, " \
                    "we recommend you start by importing your information from " \
                    "an Excel/CSV file. " \
                    "After you have your file ready, you can open it by clicking the" \
                    " \"File > Open File\" button and selecting your file. " \
                    "You can also open a file by dragging and dropping it onto the window. " \
                    "After opening the file, you should see it populate the main window like a spreadsheet.</p>" \
                    "<h3>What Do You Need In the Spreadsheet?</h3>" \
                    "<p>All you need is to specify the \"Sample name\" column with the file path. " \
                    "However, there are other optional parameters that you can specifiy. Note, capitalization is" \
                    " important, so make sure to specify caps carefully.</p>"

        html_str += array_to_html(basic_parameters, cols=["Parameter", "Required", "Description"], rows=None,
                                  colors=None, index=False, sortable=False)

        html_str += "<h3>Running UPP</h3><p>" \
                    "After you have opened your file, you can run UPP by clicking the \"Run\" button.</p> " \
                    "<p>There are two options to select to speed up the deconvolution and processing. " \
                    "The first is to use the already converted data. " \
                    "If you have already converted and averaged the data into a single spectrum," \
                    " you can select this option to skip the conversion step " \
                    "and speed up the program while you optimize the deconvolution.</p><p>" \
                    "The second option is to run the deconvolution engine. " \
                    "If you have already run the deconvolution engine on your data, " \
                    "you can select this option to skip the deconvolution step and speed up the program. " \
                    "This option is most useful if you have already deconvolved the data" \
                    " and want to adjust the peak matching or analysis.</p> "

        html_str += "<h3>Outputs</h3> <p>" \
                    "After running UPP, there are two key outputs. " \
                    "First, you will see one or more new tabs appear in the main window. " \
                    "These tabs will contain the results of the deconvolution and analysis. " \
                    "The results are saved to a \"results.xlsx\" file.</p> " \
                    "<p>Second, each deconvolution will generate an HTML reports that will be saved" \
                    " in the same directory as your data file in the _unidecfiles folder. " \
                    "You can open these reports in a web browser by clicking the \"Open All HTML Reports\" button. " \
                    "You can also open individual files by double clicking on individual cells.</p> "

        html_str += "<h3>Adjusting the Deconvolution Parameters</h3> <p>" \
                    "UPP will use the default UniDec parameters for deconvolution. " \
                    "However, you can adjust the deconvolution parameters " \
                    "by adding these optional rows in the spreadsheet: </p> "

        html_str += array_to_html(config_parameters, cols=["Parameter", "Required", "Description"], rows=None,
                                  colors=None, index=False, sortable=False)

        html_str += "<h2>Advanced Features</h2> <h3>Developing Workflows</h3><p>" \
                    "After you deconvolve your data, there are lots of things you can do with it. " \
                    "Because UPP is free and open-source, you can write in new functions and features " \
                    "that are customized for your workflow.</p> <p>For example, you could read in a column " \
                    "for \"Candidate Mass\" and search the peaks to if it is present. " \
                    "Take a look at the batch.py file on GitHub for ideas." \
                    "</p> <p>If you have new ideas for recipes, feel free to reach out for help. " \
                    "We are happy to help you develop your own recipes and workflows.</p> "

        html_str += "<h3>Workflow 1: Check Correct Masses or Pairings</h3><p>" \
                    "Here is an example recipe that checks if the correct pairing of massees" \
                    " and/or sequences is present. " \
                    "The column keyword of \"Sequence {n}\" defines a protein sequence " \
                    "were {n} is a number, such as \"Sequence 1\", or some other label. " \
                    "Each sequence cell should give the amino acid sequence of the " \
                    "protein chain or the mass of the species. </p> <p> " \
                    "Another key column is \"Correct{anything}\". UPP will look for a column with \"Correct\" in it. " \
                    "The \"Correct\" column should contain the correct pairing of the protein sequences. " \
                    "For example, if you have two protein sequences, \"Sequence 1\" and \"Sequence 2\", " \
                    "the \"Correct\" column should contain the pairing of the" \
                    " two sequences written as: \"Seq1+Seq2\". " \
                    "You can also other columns like \"Incorrect Homodimer\" as a column header " \
                    "with similar definitions (Seq2+Seq2 for example) " \
                    "and UPP will check if the incorrect pairing is present. " \
                    "You can also specify \"Ignored\" columns to ignore certain pairings. " \
                    "Note, you can also specify masses directly as correct, incorrect, or ignored. </p> <p> " \
                    "Finally, you can specify a Fixed or Variable Mod File to list potential sequence " \
                    "modifications (see more info below) and a \"Tolerance\" to specify the peak matching tolerance. " \
                    "Using all this information, the workflow will then search for the correct" \
                    " and incorrectly paired masses in the deconvolution results (with any possible modifications). " \
                    "If the correct mass/pairing is present, it will color the peak green. " \
                    "If the incorrect mass/pairing is present, it will color the peak red. " \
                    "If an ignored mass/pairng is present, it will color the peak blue. " \
                    "If no mathes are found for a given peak (unknown), it will color the peak yellow. </p> <p> " \
                    "The final results spreadsheet will contain the percentage of the signal " \
                    "that is correct, incorrect, ignored, and unknown." \
                    "It will also give the percentage of correct and incorrect " \
                    "after ignoring the unknown and ignored. " \
                    "Additional details on keywords are provided below. "

        html_str += array_to_html(recipe_w, cols=["Parameter", "Required", "Description"], rows=None,
                                  colors=None, index=False, sortable=False)

        html_str += "<h3>Workflow 2: Calculate Drug-to-Antibody Ratio</h3><p>" \
                    "This recipe will calculate the drug-to-antibody ratio for an ADC. " \
                    "The column keyword of \"Protein Mass\" defines the mass of the antibody, either with a " \
                    "mass value, an amino acid sequence, or a Seq code combination. " \
                    "The column keyword of \"Drug Mass\" defines the mass of the drug. " \
                    "The column keyword of \"Max Drug\" defines the maximum number of drug molecules to consider. " \
                    "UPP will then loop through all possible combinations of drug and antibody " \
                    "and see if they match any of the detected peaks. " \
                    "If a match is found, it will color the peak green. " \
                    "If no match is found, it will color the peak yellow. " \
                    "After extracting the height for each peak, UPP will calculate the drug-to-antibody ratio. " \
                    "Additional details on keywords are provided below. "

        html_str += array_to_html(recipe_d, cols=["Parameter", "Required", "Description"], rows=None,
                                  colors=None, index=False, sortable=False)

        html_str += "<h3>How to Build Your Own Workflow</h3><p>" \
                    "Note, this will require knowledge of Python syntax," \
                    " the Pandas library, and NumPy.</p> " \
                    "<p>You can build your own workflow by modifying the batch.py file. " \
                    "First, locate the batch.py file, and then find the run_df function. " \
                    "This function first checks the column keywords to see if each workflow should be run. " \
                    "Add your own keyword check function here.</p> " \
                    "<p>Next, the run_df function loops through each row in the DataFrame. " \
                    "For each row, it will first deconvolve the data and then extract the peaks. " \
                    "The peaks are collected in a the eng.pks structure, " \
                    "which is found in unidec/modules/peakstructure.py. " \
                    "After finding the peaks, you can add your own analysis functions where it is marked. " \
                    "The key objects are the row, which contains all the parameters from the spreadsheet, " \
                    "and the eng.pks structure, which contains the peaks.</p> " \
                    "<p>I recommend that you write a function with the row as a parameter " \
                    "and that returns a newrow back. " \
                    "You can also include the pks object as a parameter, but I usually just" \
                    " get it from the engine by default. " \
                    "The row is a pandas DataFrame object, " \
                    "and you can access the different parameters with " \
                    "row[“Column Name”], for example. The peaks can be accessed with something " \
                    "like [p.mass for p in pks.peaks]. " \
                    "You can also access values like p.height or p.integral. The self.integrate variable, " \
                    "which is set by row[“Quant Mode”], " \
                    "can be accessed to select either of these.</p>" \
                    "<p>Additional peak parameters, such as p.color or p.dscore can be found in the " \
                    "peakstructure.py code. You can edit these parameters inside your " \
                    "function to edit the color of the peaks or the label with something like: " \
                    "pks.peaks[i].color = [0, 1, 0]. See the dar_calc function in " \
                    "unidec/modules/matchtool.py for a simple example. </p>" \
                    "<p>After pulling the relevant information from the row (your spreadsheet) and the pks " \
                    "(the deconvolution results), you can do any necessary calculations and write the results" \
                    " back into the original row object, using something like row[“Result”]. " \
                    "You can name it whatever you want. Importantly, unlike the pks object, " \
                    "just modifying the row inside your function will not change it in the larger program. " \
                    "To get the modified row out, your need to return the row as the output from the function " \
                    "(which I save as newrow) and then use the set_row_merge function to merge the newrow " \
                    "into the rundf. You should just be able to copy the syntax from other workflows here.</p>" \
                    "<p>If you would like to access more features of the deconvolved data than you can find " \
                    "in the peaks, you can use the self.eng, which accesses the entire UniDec engine. " \
                    "The engine code can be found under engine.py and unidec_enginebase.py. " \
                    "The key data is stored under eng.data, which is a DataContainer object " \
                    "(see modules/unidecstructure.py). For example, self.eng.data.massdat is a numpy array with " \
                    "the 1D deconvolved mass distribution, with mass in the first position " \
                    "and intensity in the second. </p>" \
                    "<p>That’s it! The row contains all the key inputs and outputs from the spreadsheet. " \
                    "The eng.pks or eng.data contain the deconvolution data. Let me know what questions you " \
                    "have and what ideas you come up with. I would love to add additional workflows, " \
                    "so if you send me your code, I can merge it in for anyone to use.</p> " \

        html_str += "</body></html>"


        self.htmlstr = html_str

        # For Writing to File
        # Copy html_str to clipboard
        '''if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(wx.TextDataObject(html_str))
            wx.TheClipboard.Close()'''

        html.SetPage(html_str)


    def open_in_browser(self, event):
        # Open in Browser
        with open("help.html", "w") as f:
            f.write(self.htmlstr)
        webbrowser.open("help.html")


class MyFileDropTarget(wx.FileDropTarget):
    """"""

    def __init__(self, window):
        """Constructor"""
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, filenames):
        """
        When files are dropped, either open a single file or run in batch.
        """
        path = filenames[0]
        self.window.load_file(path)
        return 0


class UPPApp(wx.Frame):
    """"""

    def __init__(self, nrows=2, ncolumns=2, title="UniDec Processing Pipeline"):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title=title, size=(1800, 600))
        self.use_decon = True
        self.use_converted = False
        self.use_interactive = False
        self.make_combined_peaks = True
        self.allpng = False
        self.bpeng = BPEngine(parent=self)

        try:
            if os.path.isfile(self.bpeng.eng.config.iconfile):
                favicon = wx.Icon(self.bpeng.eng.config.iconfile, wx.BITMAP_TYPE_ANY)
                wx.Frame.SetIcon(self, favicon)
                self.icon_path = os.path.abspath(self.bpeng.eng.config.iconfile)
            else:
                self.icon_path = None
        except Exception as e:
            print(e)
            self.icon_path = None

        menu = wx.Menu()
        # Open File Menu
        open_file_menu_item = menu.Append(wx.ID_ANY, "Open File", "Open a CSV or Excel file")
        self.Bind(wx.EVT_MENU, self.on_load_file, open_file_menu_item)
        menu.AppendSeparator()

        # Save File Menu
        save_file_menu_item = menu.Append(wx.ID_ANY, "Save File", "Save a CSV or Excel file")
        self.Bind(wx.EVT_MENU, self.on_save_file, save_file_menu_item)
        menu.AppendSeparator()

        # Add Files Menu
        add_files_menu_item = menu.Append(wx.ID_ANY, "Add Data Files", "Add Data Files")
        self.Bind(wx.EVT_MENU, self.on_add_files, add_files_menu_item)

        # Clear everything on the panel
        clear_everything_menu_item = menu.Append(wx.ID_ANY, "Clear All", "Clear Everything")
        self.Bind(wx.EVT_MENU, self.clear_all, clear_everything_menu_item)

        help_menu = wx.Menu()
        # Open File Menu
        help_manu_item = help_menu.Append(wx.ID_ANY, "Help Me!", "Open a help page")
        self.Bind(wx.EVT_MENU, self.on_help_page, help_manu_item)

        # Create the menubar
        menuBar = wx.MenuBar()
        menuBar.Append(menu, "&File")
        menuBar.Append(help_menu, "&Help")
        self.SetMenuBar(menuBar)

        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Insert a button and bind it with a handler called on_run
        self.runbtn = wx.Button(panel, label="Run All")
        self.runbtn.Bind(wx.EVT_BUTTON, self.on_run)
        hsizer.Add(self.runbtn, 0)

        # Insert a button and bind it with a handler called on_run
        self.runbtn2 = wx.Button(panel, label="Run Selected")
        self.runbtn2.Bind(wx.EVT_BUTTON, self.on_run_selected)
        hsizer.Add(self.runbtn2, 0)

        # Insert Spacer Text
        hsizer.Add(wx.StaticText(panel, label="     "), 0)

        # Insert a checkbox to select whether to use already converted data
        self.useconvbox = wx.CheckBox(panel, label="Use Converted Data  ")
        hsizer.Add(self.useconvbox, 0, wx.EXPAND)
        self.useconvbox.SetValue(self.use_converted)

        # Insert a checkbox to select whether to use already deconvolved data
        self.usedeconbox = wx.CheckBox(panel, label="Deconvolve Data  ")
        hsizer.Add(self.usedeconbox, 0, wx.EXPAND)
        self.usedeconbox.SetValue(self.use_decon)

        # Insert a checkbox to select whether to generate interactive HTML reports
        self.interactivebox = wx.CheckBox(panel, label="Interactive Reports  ")
        hsizer.Add(self.interactivebox, 0, wx.EXPAND)
        self.interactivebox.SetValue(self.use_interactive)

        # Insert a checkbox to select whether to generate HTML reports with all PNG images and not SVG. Saves a bit of space.
        self.pngbox = wx.CheckBox(panel, label="PNG Figures  ")
        hsizer.Add(self.pngbox, 0, wx.EXPAND)
        self.pngbox.SetValue(self.allpng)

        # Insert a checkbox to select whether to generate a combined peak list
        self.peaklistbox = wx.CheckBox(panel, label="Generate Combined Peak List  ")
        hsizer.Add(self.peaklistbox, 0, wx.EXPAND)
        self.peaklistbox.SetValue(self.make_combined_peaks)

        # Insert Spacer Text
        hsizer.Add(wx.StaticText(panel, label="   "), 0)

        # Insert a button for Open Global HTML Reports and bind to function
        btn = wx.Button(panel, label="Open Combined Report")
        btn.Bind(wx.EVT_BUTTON, self.on_open_global_html)
        btn.SetBackgroundColour("#5B2C6F")
        btn.SetOwnForegroundColour("#FFFFFF")
        hsizer.Add(btn, 0)

        # Insert a button for Open All HTML Reports and bind to function
        btn = wx.Button(panel, label="Open All Reports")
        btn.Bind(wx.EVT_BUTTON, self.on_open_all_html)
        btn.SetBackgroundColour("#F15628")
        btn.SetOwnForegroundColour("#0B31A5")
        hsizer.Add(btn, 0)

        # Insert a button for Open the Results File and bind to function
        btn = wx.Button(panel, label="Open Results")
        btn.Bind(wx.EVT_BUTTON, self.open_results_file)
        btn.SetBackgroundColour("#196F3D")
        btn.SetOwnForegroundColour("#FFFFFF")
        hsizer.Add(btn, 0)

        # Insert a button for Run in UniDec and bind to function
        btn = wx.Button(panel, label="Open in UniDec")
        btn.Bind(wx.EVT_BUTTON, self.on_open_unidec)
        btn.SetBackgroundColour("#1F618D")
        btn.SetOwnForegroundColour("#FFFF00")
        hsizer.Add(btn, 0)

        # Insert Spacer Text
        # hsizer.Add(wx.StaticText(panel, label="   "), 0)
        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        # Insert a button to hide columns
        self.hidebtn = wx.Button(panel, label="Hide Columns")
        self.hidebtn.Bind(wx.EVT_BUTTON, self.on_hide_columns)
        hsizer2.Add(self.hidebtn, 0)
        self.hide_col_flag = False

        hsizer2.Add(wx.StaticText(panel, label="   "), 0)

        # Insert a button to hide columns with height in the title
        self.hideheightbtn = wx.Button(panel, label="Hide Height Columns")
        self.hideheightbtn.Bind(wx.EVT_BUTTON, self.on_hide_height_columns)
        self.hideheightbtn.SetBackgroundColour("#5D6D7E")
        self.hideheightbtn.SetOwnForegroundColour("#FFFFFF")
        hsizer2.Add(self.hideheightbtn, 0)
        self.hide_height_flag = False

        # Insert a button to hide columns with % in the title
        self.hidepercentbtn = wx.Button(panel, label="Hide % Columns")
        self.hidepercentbtn.Bind(wx.EVT_BUTTON, self.on_hide_percent_columns)
        self.hidepercentbtn.SetBackgroundColour("#D5D8DC")
        # self.hidepercentbtn.SetOwnForegroundColour("#FFFFFF")
        hsizer2.Add(self.hidepercentbtn, 0)
        self.hide_percentcol_flag = False

        # Insert Spacer Text
        hsizer2.Add(wx.StaticText(panel, label="   "), 0)

        # Insert a button to hide columns that are empty
        self.hideemptybtn = wx.Button(panel, label="Hide Empty Columns")
        self.hideemptybtn.Bind(wx.EVT_BUTTON, self.on_hide_empty_columns)
        hsizer2.Add(self.hideemptybtn, 0)

        sizer.Add(hsizer, 0, wx.ALL | wx.EXPAND)

        self.ss = SpreadsheetPanel(self, panel, nrows, ncolumns).ss
        self.ss.set_col_headers(["Sample name", "Data Directory"])
        sizer.Add(self.ss, 1, wx.EXPAND)

        sizer.Add(hsizer2, 0, wx.ALL | wx.EXPAND)

        file_drop_target = MyFileDropTarget(self)
        self.ss.SetDropTarget(file_drop_target)
        # self.ss.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.on_cell_clicked)

        panel.SetSizer(sizer)

        # Create Status bar
        self.CreateStatusBar(3)

        self.Show()

    def on_run(self, event=None):
        self.SetStatusText("Running...", 2)
        self.runbtn.SetBackgroundColour("red")
        self.get_from_gui()
        wx.Yield()
        self.bpeng.run_df(decon=self.use_decon, use_converted=self.use_converted, interactive=self.use_interactive,
                          write_peaks=self.make_combined_peaks, allpng=self.allpng)
        self.load_to_gui()
        self.runbtn.SetBackgroundColour("green")
        if not self.hide_col_flag:
            self.on_hide_columns()

        self.SetStatusText("Outputs: " + self.bpeng.outfile, 1)
        self.SetStatusText("Completed: " + str(self.bpeng.runtime) + " seconds", 2)

    def on_run_selected(self, event=None, rows=None):
        # Get from GUI
        self.SetStatusText("Running...", 2)
        self.runbtn2.SetBackgroundColour("red")
        if rows is None:
            # Get Selected Rows
            selected_rows = list(self.ss.get_selected_rows())
        else:
            selected_rows = rows
        print("Running Selected Rows:", selected_rows)
        # Get Sub Dataframe with Selected Rows
        self.get_from_gui()
        wx.Yield()

        # RUN THE SELECTED ROWS
        topdf = deepcopy(self.bpeng.rundf)
        subdf = self.bpeng.rundf.iloc[selected_rows]
        toppeaks = deepcopy(self.bpeng.pks)
        # Run SubDF
        subdf2 = self.bpeng.run_df(df=subdf, decon=self.use_decon, use_converted=self.use_converted,
                                   interactive=self.use_interactive, write_xlsx=False, write_html=False,
                                   write_peaks=False, allpng=self.allpng)

        # Update the main dataframe
        # topdf.iloc[selected_rows] = subdf2
        topdf = set_row_merge(topdf, subdf, selected_rows)
        self.bpeng.rundf = topdf

        # Write the results
        self.bpeng.write_xlsx()

        if self.make_combined_peaks:
            if toppeaks is not None:
                # Merge in the new peaks
                toppeaks = toppeaks.merge_in_peaks(self.bpeng.pks)
                # Write the peaks
                self.bpeng.pks = toppeaks
                self.bpeng.write_peaks()
            else:
                print("Missing peaks for the whole data set. Outputting only the selected peaks.")
                self.bpeng.write_peaks()

        # Load to GUI
        self.load_to_gui()
        # Finish by coloring the button green
        self.runbtn2.SetBackgroundColour("green")
        if not self.hide_col_flag:
            self.on_hide_columns()

        self.SetStatusText("Outputs: " + self.bpeng.outfile, 1)
        self.SetStatusText("Completed: " + str(self.bpeng.runtime) + " seconds", 2)

    def clear_all(self, event=None):
        self.ss.delete_all()
        self.ss.set_col_headers(["Sample name", "Data Directory"])

    def load_file(self, filename):
        print("Loading File:", filename)
        # Set status
        self.SetStatusText("Inputs: " + filename, 0)
        try:
            self.ss.delete_all()
        except Exception:
            pass
        self.bpeng.top_dir = os.path.dirname(filename)
        self.bpeng.filename = filename
        self.bpeng.rundf = file_to_df(filename)
        self.load_to_gui()
        # dirname = os.path.dirname(filename)
        # self.set_dir_tet_box(dirname)
        self.reset_hidden_columns()

    def on_load_file(self, event):
        print("Load button pressed")
        # Create a file dialog
        with wx.FileDialog(self, "Open CSV or Excel File",
                           wildcard="CSV or Excel files (*.csv; *.xlsx; *.xls)|*.csv; *.xlsx; *.xls|"
                                    "CSV files (*.csv)|*.csv|"
                                    "Excel files (*.xlsx; *.xls)|*.xlsx; *.xls",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            # Show the dialog and retrieve the user response. If it is the OK response,
            # process the data.
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            self.load_file(pathname)

    def on_save_file(self, event):
        print("Save button pressed")
        # Create a file dialog
        with wx.FileDialog(self, "Save CSV or Excel File",
                           wildcard="CSV or Excel files (*.csv; *.xlsx; *.xls)|*.csv; *.xlsx; *.xls|"
                                    "CSV files (*.csv)|*.csv|"
                                    "Excel files (*.xlsx; *.xls)|*.xlsx; *.xls",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
            # Show the dialog and retrieve the user response. If it is the OK response,
            # process the data.
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            self.ss.save_file(pathname)

    # def set_dir_tet_box(self, dirname):
    #    self.dirtxtbox.SetValue(dirname)

    def get_from_gui(self):
        self.use_converted = self.useconvbox.GetValue()
        self.use_decon = self.usedeconbox.GetValue()
        self.use_interactive = self.interactivebox.GetValue()
        self.make_combined_peaks = self.peaklistbox.GetValue()
        self.allpng = self.pngbox.GetValue()

        self.ss.remove_empty()
        ssdf = self.ss.get_df()
        self.bpeng.rundf = ssdf

    def load_to_gui(self):
        self.ss.set_df(self.bpeng.rundf)
        self.color_columns()

    def open_results_file(self, event):
        print("Open Results button pressed")
        # Open Results File in Explorer
        if os.path.isfile(self.bpeng.outfile):
            os.startfile(self.bpeng.outfile)

    def on_open_all_html(self, event):
        print("Open All HTML Reports button pressed")
        self.bpeng.open_all_html()

    def on_open_global_html(self, event):
        print("Open Global HTML Report button pressed")
        self.bpeng.open_global_html()

    def on_open_unidec(self, event):
        ssdf = self.ss.get_df()
        self.bpeng.rundf = ssdf
        selected_rows = list(self.ss.get_selected_rows())
        print(selected_rows)
        row = self.bpeng.rundf.iloc[selected_rows[0]]
        self.open_unidec(row)

    def open_unidec(self, row):
        print("Opening in UniDec:", row)
        filepath = self.bpeng.get_file_path(row)
        if filepath is not None:
            print("Launching UniDec:", filepath)
            app = UniDecApp(path=filepath)
            app.eng.unidec_imports(efficiency=False)
            app.after_unidec_run()
            app.on_pick_peaks()
            if self.bpeng.correct_pair_mode:
                self.bpeng.run_correct_pair(row, app.eng.pks)
                app.after_pick_peaks()
            elif self.bpeng.dar_mode:
                self.bpeng.run_dar(row, app.eng.pks)
                app.after_pick_peaks()
            app.start()

    def on_add_files(self, event=None):

        wildcard = "CSV or Excel files (*.csv; *.xlsx; *.xls)|*.csv; *.xlsx; *.xls|CSV files (*.csv)|*.csv|Excel files (*.xlsx; *.xls)|*.xlsx; *.xls"
        wildcard = "Any files (*.*) |*.*| " \
                   "Known file types (*.raw; *.d; *.mzML; *.mzXML; *.txt; *.csv; *.dat; *.npz)|" \
                   "*.raw; *.d; *.mzML; *.mzXML; *.txt; *.csv; *.dat; *.npz|" \
                   "Thermo RAW files (*.raw)|*.raw|" \
                   "Agilent D files (*.d)|*.d|" \
                   "mzML files (*.mzML)|*.mzML|" \
                   "mzXML files (*.mzXML)|*.mzXML|" \
                   "Text files (*.txt)|*.txt|" \
                   "CSV files (*.csv)|*.csv|" \
                   "Dat files (*.dat)|*.dat|" \
                   "NPZ files (*.npz)|*.npz"

        # Create a file selection dialog
        with wx.FileDialog(self, "Select Files to Add",
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE) as fileDialog:
            # Show the dialog and retrieve the user response. If it is the OK response,
            # process the data.

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            # Proceed loading the file chosen by the user
            paths = fileDialog.GetPaths()
            self.add_files(paths)

    def add_files(self, paths):
        print("Adding Files:", paths)
        self.get_from_gui()
        # Add paths list to the datafram in the "Sample name" column
        sample_names = [os.path.basename(path) for path in paths]
        data_dir = [os.path.dirname(path) for path in paths]
        newdf = pd.DataFrame({"Sample name": sample_names, "Data Directory": data_dir})
        self.bpeng.rundf = pd.concat([self.bpeng.rundf, newdf], ignore_index=True)

        self.load_to_gui()
        self.reset_hidden_columns()

    def on_hide_columns(self, event=None, reset=False):
        columns_to_hide = ["Tolerance", "File", "Time", "Config", "Sequence", "Directory", "Matches",
                           "Global Fixed Mod", "Favored Match", "Apply Fixed Mods", "Disulfides Oxidized"]
        if not self.hide_col_flag and not reset:
            for keyword in columns_to_hide:
                self.ss.hide_columns_by_keyword(keyword)
            self.hide_col_flag = True
            self.hidebtn.SetLabel("Show Columns")
        else:
            self.hidebtn.SetLabel("Hide Columns")
            self.ss.show_all_columns()
            self.hide_col_flag = False

    def reset_hidden_columns(self):
        self.on_hide_columns(reset=True)
        self.on_hide_height_columns(reset=True)
        self.on_hide_percent_columns(reset=True)

    def on_hide_empty_columns(self, event=None):
        self.ss.hide_empty_columns()

    def on_hide_height_columns(self, event=None, reset=False):
        if not self.hide_height_flag and not reset:
            self.ss.hide_columns_by_keyword("Height")
            self.hide_height_flag = True
            self.hideheightbtn.SetLabel("Show Height Columns")
        else:
            self.hideheightbtn.SetLabel("Hide Height Columns")
            self.ss.show_columns_by_keyword("Height")
            self.hide_height_flag = False

    def on_hide_percent_columns(self, event=None, reset=False):
        if not self.hide_percentcol_flag and not reset:
            self.ss.hide_columns_by_keyword(" %")
            self.hide_percentcol_flag = True
            self.hidepercentbtn.SetLabel("Show % Columns")
        else:
            self.hidepercentbtn.SetLabel("Hide % Columns")
            self.ss.show_columns_by_keyword(" %")
            self.hide_percentcol_flag = False

    def color_columns(self):
        for row in basic_parameters:
            if row[1]:
                self.ss.color_columns_by_keyword(row[0], "#D1F2EB")
            else:
                self.ss.color_columns_by_keyword(row[0], "#D6EAF8")

        for row in config_parameters:
            if row[1]:
                self.ss.color_columns_by_keyword(row[0], "#FEF9E7")
            else:
                self.ss.color_columns_by_keyword(row[0], "#FEF9E7")

        for row in recipe_w:
            if row[1]:
                self.ss.color_columns_by_keyword(row[0], "#BB8FCE")
            else:
                self.ss.color_columns_by_keyword(row[0], "#E8DAEF")

        for row in recipe_d:
            if row[1]:
                self.ss.color_columns_by_keyword(row[0], "#D98880")
            else:
                self.ss.color_columns_by_keyword(row[0], "#F5B7B1")

        self.ss.color_columns_by_keyword("Height", "#5D6D7E")
        self.ss.color_columns_by_keyword("%", "#D5D8DC")
        self.ss.color_columns_by_keyword("BsAb Pairing Calculated", "#7DCEA0")
        self.ss.color_columns_by_keyword("DAR", "#7DCEA0")
        self.ss.color_columns_by_keyword("Light Chain Scrambled", "#F7DC6F")

        self.ss.color_columns_by_keyword("Sequence", "#FDEBD0")
        self.ss.color_columns_by_keyword("Matches", "#D5F5E3")

    def update_progress(self, i, total_n):
        status_text = "Processing file {} of {}".format(i + 1, total_n)
        self.SetStatusText(status_text, 2)
        # update gui
        wx.Yield()

    """
    def on_cell_clicked(self, event=None):
        print("Cell clicked")"""

    def on_help_page(self, event=None):
        print("Help button pressed")
        dlg = HelpDlg()
        dlg.Show()

    def on_exit(self, event=None):
        self.Close()


if __name__ == "__main__":
    app = wx.App()
    frame = UPPApp()
    #frame.usedeconbox.SetValue(False)
    path = "C:\\Data\\Wilson_Genentech\\sequences_short.xlsx"
    path = "C:\\Data\\Wilson_Genentech\\BsAb\\BsAb test short.xlsx"
    path = "C:\\Data\\UPPDemo\\BsAb\\BsAb test - Copy.xlsx"
    # path = "C:\\Data\\UPPDemo\\DAR\\Biotin UPP template WP_MTM_DoubleDec.xlsx"
    # path = "C:\\Data\\Wilson_Genentech\\BsAb\\BsAb test.xlsx"
    # path = "C:\\Data\\Wilson_Genentech\\DAR\\Biotin UPP template test.xlsx"
    # path = "C:\\Data\\UPPDemo\\BsAb\\outliers.xlsx"
    # path = "C:\\Data\\UPPDemo\\BsAb\\BsAb test short.xlsx"
    path = "C:\\Data\\UPPDemo\\DAR\\Biotin UPP template WP_MTM.xlsx"
    path = "C:\\Data\\UPPDemo\\FileOpening\\filetest.xlsx"
    path = "C:\\Data\\Wilson_Genentech\\BsAb\\BsAb test short_repeat.xlsx"
    # frame.on_help_page()
    # exit()
    if False and platform.node() == 'MTM-VOSTRO':
        frame.useconvbox.SetValue(False)
        frame.usedeconbox.SetValue(True)
        frame.load_file(path)
        # frame.set_dir_tet_box("C:\\Data\\Wilson_Genentech\\Data")
        # print(df)
        frame.on_run()
        # frame.on_run_selected(rows=[1])
        # frame.on_run_selected(rows=[0])
        # frame.on_add_files()

    app.MainLoop()
