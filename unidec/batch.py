from unidec.engine import UniDec
from unidec.modules.matchtools import *
from unidec.tools import known_extensions
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
    return None


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
        print("Data Directory:", self.data_dir)
        starttime = time.perf_counter()
        if df is not None:
            self.rundf = df

        percs = []
        percs2 = []
        matches = []
        htmlfiles = []

        self.correct_pair_mode = check_for_correct_in_keys(self.rundf)

        for i, row in self.rundf.iterrows():
            file = row["Sample name"]
            if "Tolerance (Da)" in row:
                try:
                    self.tolerance = float(row["Tolerance (Da)"])
                except Exception:
                    pass
            if "Mod File" in row:
                self.modfile = row["Mod File"]
                if os.path.isfile(self.modfile):
                    self.moddf = file_to_df(self.modfile)

            if "Start Time" in row and "End Time" in row:
                time_range = [row["Start Time"], row["End Time"]]
            else:
                time_range = None

            if "Data Directory" in row:
                self.data_dir = row["Data Directory"]

            path = find_file(file, self.data_dir, use_converted)
            if os.path.isfile(path):
                print("Opening:", path)
                self.eng.open_file(path, time_range=time_range)
                #self.eng.config.peakwindow = 15
                #self.eng.config.peakthresh = 0.01
                #self.eng.config.massub = 200000

                self.eng = set_param_from_row(self.eng, row)

                if decon:
                    self.eng.autorun()
                else:
                    self.eng.unidec_imports(efficiency=True)
                    self.eng.pick_peaks()

                if self.correct_pair_mode:
                    perc, perc2, matcharray, matchstring = UPP_check_peaks(self.rundf.iloc[i], self.eng.pks,
                                                                           self.tolerance,self.moddf)
                    percs.append(perc)
                    percs2.append(perc2)
                    matches.append(matchstring)

                outfile = self.eng.gen_html_report(open_in_browser=False)
                htmlfiles.append(outfile)
            else:
                print("File not found:", path)
                if self.correct_pair_mode:
                    percs.append([0, 0, 0])
                    percs2.append([0, 0])
                    matches.append("")
                htmlfiles.append("")

        if self.correct_pair_mode:
            percs = np.array(percs)
            percs2 = np.array(percs2)

            self.rundf["Correct %"] = percs[:, 0]
            self.rundf["Incorrect %"] = percs[:, 1]
            self.rundf["Unmatched %"] = percs[:, 2]
            self.rundf["Correct % Matched Only"] = percs2[:, 0]
            self.rundf["Incorrect % Matched Only"] = percs2[:, 1]
            self.rundf["Matches"] = matches
        self.rundf["Reports"] = htmlfiles

        if self.top_dir == "" and self.data_dir != "":
            self.top_dir = os.path.dirname(self.data_dir)

        self.rundf.to_excel(os.path.join(self.top_dir, "results.xlsx"))
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
