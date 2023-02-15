from unidec.engine import UniDec
from unidec.modules.matchtools import *
from unidec.tools import known_extensions, strip_char_from_string
import os
import numpy as np
import time
import webbrowser


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
    return eng


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
                     ["Config Low Mass", False, "Deconvolution Setting: The Low Mass Limit in Da"]]

recipe_w = [["Tolerance (Da)", False, "The Tolerance in Da. Default is 50 Da if not specified."],
            ["Mod File", False, "The File Name or Path of the Mod File. "
                                "Can be either Excel or CSV. The file should have \"Mass\" and \"Name\" "
                                "columns with these exact headers. If not specified, no modifications will be used."],
            ["Sequence {n}", True, "The amino acid sequence for the {n}th protein. Can be multiple sequences, "
                                   "each as their own columns."],
            ["Correct", True, "The Correct Pairing. This should be a list of the correct pairs, "
                               "listed as Seq1+Seq2, for example."],
            ["{any other keywords}", True, "Incorrect Pairings. This should be a list of the incorrect pairs, "
                                            "listed as Seq2+Seq2, for example. "
                                            "The column header can be anything as long as \"Seq\" is in the cell."]
            ]


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

    def run_file(self, file=None, decon=True, use_converted=True):
        self.top_dir = os.path.dirname(file)
        self.rundf = file_to_df(file)
        self.run_df(decon=decon, use_converted=use_converted)

    def run_df(self, df=None, decon=True, use_converted=True):
        # Print the data directory and start the clock
        print("Data Directory:", self.data_dir)
        starttime = time.perf_counter()
        # Set the Pandas DataFrame
        if df is not None:
            self.rundf = df

        # Set up the lists for the results
        percs = []
        percs2 = []
        matches = []
        htmlfiles = []

        # Check if the "Correct" column is in the DataFrame to start correct pair mode
        self.correct_pair_mode = check_for_correct_in_keys(self.rundf)

        # Loop through the DataFrame
        for i, row in self.rundf.iterrows():
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

            # Look for the time range
            starttime = 0
            endtime = 1e12
            if "Start Time" in row:
                starttime = row["Start Time"]
            if "End Time" in row:
                endtime = row["End Time"]
            if starttime != 0 or endtime != 1e12:
                time_range = [startime, endtime]
            else:
                time_range = None

            # Find the file
            path = find_file(file, self.data_dir, use_converted)

            # If the file exists, open it
            if os.path.isfile(path):
                print("Opening:", path)
                self.eng.open_file(path, time_range=time_range)

                # Set the deconvolution parameters from the DataFrame
                self.eng = set_param_from_row(self.eng, row)

                # Run the deconvolution or import the prior deconvolution results
                if decon:
                    self.eng.autorun()
                else:
                    self.eng.unidec_imports(efficiency=True)
                    self.eng.pick_peaks()

                # The First Recipe, correct pair mode
                if self.correct_pair_mode:
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
                    perc, perc2, matcharray, matchstring = UPP_check_peaks(self.rundf.iloc[i], self.eng.pks,
                                                                           self.tolerance, self.moddf)
                    # Fill in the results
                    percs.append(perc)
                    percs2.append(perc2)
                    matches.append(matchstring)

                # Generate the HTML report
                outfile = self.eng.gen_html_report(open_in_browser=False)
                htmlfiles.append(outfile)
            else:
                # When files are not found, print the error and add empty results
                print("File not found:", path)
                if self.correct_pair_mode:
                    percs.append([0, 0, 0])
                    percs2.append([0, 0])
                    matches.append("")
                htmlfiles.append("")

        # Set the results array in the DataFrame
        if self.correct_pair_mode:
            percs = np.array(percs)
            percs2 = np.array(percs2)

            self.rundf["Correct %"] = percs[:, 0]
            self.rundf["Incorrect %"] = percs[:, 1]
            self.rundf["Unmatched %"] = percs[:, 2]
            self.rundf["Correct % Matched Only"] = percs2[:, 0]
            self.rundf["Incorrect % Matched Only"] = percs2[:, 1]
            self.rundf["Matches"] = matches

        # Set the HTML file names in the DataFrame
        self.rundf["Reports"] = htmlfiles

        # If the top directory is not set, set it as the directory above the data directory
        if self.top_dir == "" and self.data_dir != "":
            self.top_dir = os.path.dirname(self.data_dir)

        # Write the results to an Excel file to the top directory
        self.rundf.to_excel(os.path.join(self.top_dir, "results.xlsx"))

        # Print Run Time
        print("Batch Run Time:", time.perf_counter() - starttime)

    def open_all_html(self):
        for i, row in self.rundf.iterrows():
            if os.path.isfile(row["Reports"]):
                webbrowser.open(row["Reports"])


if __name__ == "__main__":
    # Example usage
    # Create a dataframe with the following columns:
    # Sample name, Mass, Charge, Tolerance (Da), Mod File, Start Time, End Time, Data Directory
    # Then run:
    # batch = UniDecBatchProcessor()
    # batch.run_df(df)
    # batch.open_all_html()
    batch = UniDecBatchProcessor()
    path = "C:\\Data\\Wilson_Genentech\\sequences_short.xlsx"
    batch.run_file(path, decon=True, use_converted=True)
    batch.open_all_html()
