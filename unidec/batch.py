from unidec.engine import UniDec
from unidec.modules.matchtools import *
from unidec.tools import known_extensions, strip_char_from_string
import os
import numpy as np
import time
import webbrowser
import sys

basic_parameters = [["Sample name", True, "The File Name or Path. File extensions are optional."],
                    ["Data Directory", False, "The directory of the data files. If you do not specify this, "
                                              "you will need to put a full path in the \"Sample name\"."],
                    ["Start Time", False, "Deconvolution Setting: The Start Time to Start Extracting "
                                          "From if Chromatography Data is Used. If not specified, "
                                          "will assume start with scan 1."],
                    ["End Time", False, "Deconvolution Setting: The End Time to Stop Extracting "
                                        "From if Chromatography Data is Used. If not specified, "
                                        "will continue to the last scan. It will sum from Start Time to End Time"]
                    ]

config_parameters = [["Config Peak Thres", False, "Deconvolution Setting: The Peak Detection Threshold"],
                     ["Config Peak Window", False, "Deconvolution Setting: The Peak Detection Window in Da"],
                     ["Config High Mass", False, "Deconvolution Setting: The High Mass Limit in Da"],
                     ["Config Low Mass", False, "Deconvolution Setting: The Low Mass Limit in Da"],
                     ["Config High m/z", False, "Deconvolution Setting: The High m/z Limit"],
                     ["Config Low m/z", False, "Deconvolution Setting: The Low m/z Limit"]
                     ]

recipe_w = [["Tolerance (Da)", False, "The Tolerance in Da. Default is 50 Da if not specified."],
            ["Mod File", False, "The File Name or Path of the Mod File. "
                                "Can be either Excel or CSV. The file should have \"Mass\" and \"Name\" "
                                "columns with these exact headers. If not specified, no modifications will be used."],
            ["Favored Match", False, "If there are several possible matches within tolerance, which to select. "
                                     "Default is \"Closest\" for the closest absolute mass. Also accepts \"Incorrect\" "
                                     "to select an incorrect mass within the tolerance even if a correct mass is "
                                     "closesr. Other keywords like \"Ignore\" and \"Correct\" can be used to select "
                                     "specific types of matches over others."],
            ["Sequence {n}", True, "The amino acid sequence or mass for the {n}th protein. Can be multiple sequences, "
                                   "each as their own columns. If it can convert to float, it will. "
                                   "If not, it will assume it is an amino acid sequence with 1 letter codes."],
            ["Correct", True, "The Correct Pairing. This should be a list of the correct pairs, "
                              "listed as Seq1+Seq2, for example."],
            ["Ignore", False, "Known species that can be ignored. This should be similar to the list of correct pairs, "
                              "listed as Seq1+Seq2, for example. Can be Seq5 if you want to just ignore Seq5. "
                              "Note, mods will also be applied to this."],
            ["Incorrect", True, "Incorrect Pairings. This should be a list of the incorrect pairs, "
                                           "listed as Seq2+Seq2, for example. "]
            ]


def find_file(fname, folder, use_converted=True):
    if use_converted:
        extensions = known_extensions[::-1]
    else:
        extensions = known_extensions
    for ext in extensions:
        if os.path.exists(os.path.join(folder, fname + ext)):
            return os.path.join(folder, fname + ext)
    return fname


def check_for_correct_in_keys(df):
    for k in df.keys():
        if "Correct" in k:
            return True
    return False


def set_param_from_row(eng, row):
    for k in row.keys():
        if "Config Peak Window" in k:
            try:
                eng.config.peakwindow = float(row[k])
            except Exception as e:
                print("Error setting peak window", k, row[k], e)
        if "Config Peak Thres" in k:
            try:
                eng.config.peakthresh = float(row[k])
            except Exception as e:
                print("Error setting peak threshold", k, row[k], e)
        if "Config Low Mass" in k:
            try:
                eng.config.masslb = float(row[k])
            except Exception as e:
                print("Error setting low mass", k, row[k], e)
        if "Config High Mass" in k:
            try:
                eng.config.massub = float(row[k])
            except Exception as e:
                print("Error setting high mass", k, row[k], e)
        if "Config Low m/z" in k or "Config Low mz" in k:
            try:
                eng.config.minmz = float(row[k])
            except Exception as e:
                print("Error setting low m/z", k, row[k], e)
        if "Config High m/z" in k or "Config High mz" in k:
            try:
                eng.config.maxmz = float(row[k])
            except Exception as e:
                print("Error setting high m/z", k, row[k], e)

        # print(eng.config.maxmz, eng.config.minmz, k)
    return eng


def get_time_range(row):
    # Look for the time range
    starttime = 0
    endtime = 1e12
    if "Start Time" in row:
        starttime = row["Start Time"]
    if "End Time" in row:
        endtime = row["End Time"]
    if starttime != 0 or endtime != 1e12:
        try:
            starttime = float(starttime)
            endtime = float(endtime)
            time_range = [starttime, endtime]
        except Exception as e:
            print("Error setting time range", starttime, endtime, e)
            time_range = None
    else:
        time_range = None
    return time_range


def set_row_merge(topdf, subdf, indexes):
    if subdf.ndim == 1:
        subdf = subdf.to_frame().transpose()

    for index in indexes:
        for k in subdf.keys():
            if k in topdf.keys():
                topdf.loc[index, k] = subdf[k][index]
                pass
            else:
                topdf[k] = subdf[k]
    return topdf


def remove_columns(df, key):
    for k in df.keys():
        if key in k:
            df = df.drop(k, axis=1)
            print("Dropping", k)
    return df


class UniDecBatchProcessor(object):
    def __init__(self):
        self.eng = UniDec()
        self.tolerance = 50
        self.data_dir = ""
        self.top_dir = ""
        self.rundf = None
        self.modfile = None
        self.moddf = None
        self.correct_pair_mode = False
        self.time_range = None

    def run_file(self, file=None, decon=True, use_converted=True, interactive=False):
        self.top_dir = os.path.dirname(file)
        self.rundf = file_to_df(file)
        self.run_df(decon=decon, use_converted=use_converted, interactive=interactive)

    def run_df(self, df=None, decon=True, use_converted=True, interactive=False):

        # Print the data directory and start the clock
        print("Data Directory:", self.data_dir)
        clockstart = time.perf_counter()
        # Set the Pandas DataFrame
        if df is not None:
            self.rundf = df

        # Set up the lists for the results
        htmlfiles = []

        # Check if the "Correct" column is in the DataFrame to start correct pair mode
        self.correct_pair_mode = check_for_correct_in_keys(self.rundf)
        #if self.correct_pair_mode:
        #    self.rundf = remove_columns(self.rundf, "Height")

        # Loop through the DataFrame
        for i, row in self.rundf.iterrows():
            path = self.get_file_path(row, use_converted=use_converted)

            # Get the time range
            self.time_range = get_time_range(row)

            # If the file exists, open it
            if os.path.isfile(path):
                print("Opening:", path)
                self.eng.open_file(path, time_range=self.time_range)

                # Set the deconvolution parameters from the DataFrame
                self.eng = set_param_from_row(self.eng, row)

                # Run the deconvolution or import the prior deconvolution results
                if decon:
                    self.eng.autorun()
                else:
                    self.eng.unidec_imports(efficiency=False)
                    self.eng.pick_peaks()

                # The First Recipe, correct pair mode
                if self.correct_pair_mode:
                    # Run correct pair mode
                    newrow = self.run_correct_pair(row)

                    # Merge the row back in the df
                    self.rundf = set_row_merge(self.rundf, newrow, [i])

                # Generate the HTML report
                outfile = self.eng.gen_html_report(open_in_browser=False, interactive=interactive)
                htmlfiles.append(outfile)
            else:
                # When files are not found, print the error and add empty results
                print("File not found:", path)
                htmlfiles.append("")

        # Set the HTML file names in the DataFrame
        self.rundf["Reports"] = htmlfiles

        # If the top directory is not set, set it as the directory above the data directory
        if self.top_dir == "" and self.data_dir != "":
            self.top_dir = os.path.dirname(self.data_dir)

        # Write the results to an Excel file to the top directory
        outfile = os.path.join(self.top_dir, "results.xlsx")
        self.rundf.to_excel(outfile)
        print("Write to: ", outfile)
        #print(self.rundf)
        # Print Run Time
        print("Batch Run Time:", time.perf_counter() - clockstart)
        return self.rundf

    def open_all_html(self):
        for i, row in self.rundf.iterrows():
            if os.path.isfile(row["Reports"]):
                webbrowser.open(row["Reports"])

    def get_file_path(self, row, use_converted=True):
        # Read the file name
        file = row["Sample name"]
        try:
            file = strip_char_from_string(file, "\"")
        except Exception as e:
            print("Error stripping file name", file, e)

        # Look for the data directory
        if "Data Directory" in row:
            data_dir = row["Data Directory"]
            if os.path.isdir(data_dir):
                self.data_dir = data_dir
        # Find the file
        outpath = find_file(file, self.data_dir, use_converted)

        return outpath

    def run_correct_pair(self, row, pks=None):
        if pks is None:
            pks = self.eng.pks

        # Extract the tolerance for peak matching
        if "Tolerance (Da)" in row:
            try:
                self.tolerance = float(row["Tolerance (Da)"])
            except Exception:
                pass

        # Get and read the mod file
        if "Mod File" in row:
            self.modfile = row["Mod File"]
            if os.path.isfile(self.modfile):
                self.moddf = file_to_df(self.modfile)

        # Match to the correct peaks
        newrow = UPP_check_peaks(row, pks, self.tolerance, self.moddf)

        return newrow


if __name__ == "__main__":
    # Example usage
    # Create a dataframe with the following columns:
    # Sample name, Mass, Charge, Tolerance (Da), Mod File, Start Time, End Time, Data Directory
    # Then run:
    # batch = UniDecBatchProcessor()
    # batch.run_df(df)
    # batch.open_all_html()
    batch = UniDecBatchProcessor()
    if len(sys.argv) > 1:
        batch.run_file(sys.argv[1], decon=True, use_converted=True)
        batch.open_all_html()
    else:
        path = "C:\\Data\\Wilson_Genentech\\sequences_short2.xlsx"
        pd.set_option('display.max_columns', None)
        batch.run_file(path, decon=False, use_converted=True, interactive=True)

        #batch.open_all_html()
        pass
