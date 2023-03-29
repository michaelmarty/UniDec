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
                     ["Config Low m/z", False, "Deconvolution Setting: The Low m/z Limit"],
                     ["Config Sample Mass Every", False, "Deconvolution Setting: The Mass Bin Size in Da"],
                     ["Config m/z Peak FWHM", False,
                      "Deconvolution Setting: The Peak Width in m/z used for Deconvolution"],
                     ["Config m/z Peak Shape", False, "Deconvolution Setting: The peak shape function used for "
                                                      "Deconvolution. 0=gaussian, 1=lorentzian, 2=split G/L"],
                     ["Config File", False, "Path to Config File. Will load this file and use the settings. Note,"
                                            " any settings specified in the batch file will override the config file."],
                     ]

recipe_w = [["Tolerance (Da)", False, "The Tolerance in Da. Default is 50 Da if not specified."],
            ["Variable Mod File", False, "The File Name or Path of the Variable Mod File. "
                                         "Can be either Excel or CSV. The file should have \"Mass\" and \"Name\" "
                                         "columns with these exact headers. If not specified,"
                                         " no modifications will be used. "
                                         "Note, if a Variable Mod File is provided, it will apply at least one mod "
                                         "to each complex. In other words, if you only have one line in the file,  "
                                         "it will act as a global fixed mod (see below) to the whole complex. "
                                         "If you want an option for unmodified, you need to include a line with "
                                         "a mass of 0.0 and a name of \"Unmodified\" or \" \". "
             ],
            ["Fixed Mod File", False, "The File Name or Path of the Fxied Mod File. "
                                      "Can be either Excel or CSV. The file should have \"Mass\" and \"Name\" "
                                      "columns with these exact headers. Can also have \"Number\" if a "
                                      "multiple is used If not specified, no modifications will be used. "
                                      "Note, it will apply one set of all fixed mods to each sequence. "
                                      "If you specify, \"Seq1+Seq2\", it will apply the fixed mods to both sequences."],
            ["Apply Fixed Mods", False,
             "A column specifying which sequences should get the fixed mods. "
             "Should have the format of \"Seq1 Seq2 Seq3\" where Seq1 is the Sequence 1 "
             "column name, etc. Delimiters do not matter. "
             "You can specify \"All\" to apply mods to all sequences "
             "or \"None\" to apply to none. "
             "It will assume yes to all if this column is not present. "],
            ["Global Fixed Mod", False, "A column specifying a global fixed mass shift to apply to all complexes. "
                                        "Unlike Fixed Mod File, which applies a fixed modification to each sequence, "
                                        "this is applied only once to each complex. "
                                        "Also, it is a single float value rather than a file. "],
            ["Disulfides Oxidized", False,
             "A column specifying the sequences that should be fully disulfide oxidized. "
             "Should have the format of \"Seq1 Seq2 Seq3\" where Seq1 is the Sequence 1 "
             "column name, etc. Delimiters do not matter. "
             "It will assume no to all if this column is not present. "
             "You can specify \"All\" to oxidize all sequences "
             "or \"None\" to oxidize none. "
             "Will only work if sequences are amino acid codes with C. "
             "It will subtract one H mass for each C."],
            ["Favored Match", False,
             "If there are several possible matches within tolerance, which to select. "
             "Default is \"Closest\" for the closest absolute mass. Also accepts \"Incorrect\" "
             "to select an incorrect mass within the tolerance even if a correct mass is "
             "closesr. Other keywords like \"Ignore\" and \"Correct\" can be used to select "
             "specific types of matches over others."],
            ["Sequence {n}", False,
             "The amino acid sequence or mass for the {n} protein or molecule. There can be multiple sequences, "
             "each as their own columns and with unique names. If it can convert to float, it will. "
             "If not, it will assume it is an amino acid sequence with 1 letter codes."],
            ["Correct", True, "The Correct Pairing. This may be a list of the correct pairs, "
                              "listed as Seq1+Seq2, for example. Can also be a single mass value. "],
            ["Ignore", False,
             "Known species that can be ignored. This should be similar to the list of correct pairs, "
             "listed as Seq1+Seq2, for example. Can be Seq5 if you want to just ignore Seq5. "
             "Note, mods will also be applied to this. Can be a single mass value as well. "],
            ["Incorrect", False, "Incorrect Pairings. This can be a list of the incorrect pairs, "
                                 "listed as Seq2+Seq2, for example. It can also be a single mass value. "],
            ["BsAb (Correct) | LC1 Mispair (Incorrect) | LC2 Mispair (Incorrect)", False,
             "This is a special mode for bispecific antibodies. If you specify all of these three columns, "
             "it will calculate the corrected percentages for bispecific antibodies and light chain scrambling. "
             "Form more details, see https://dx.doi.org/10.1080/19420862.2016.1232217. "]
            ]

recipe_d = [
    ["Protein Mass", True, "The Protein Mass in either (1) a float in Da, (2) as an amino acid sequence, or "
                           "(3) as a string of Seq{n}+Seq{m}, where Seq{n} and Seq{m} are the names of "
                           "the columns containing the amino acid sequences or masses of individual species. "
                           "See below for more details on the Sequence {n} column. "
     ],
    ["Drug Mass", True, "The Drug Mass in Da."],
    ["Min Drug", False,
     "The minimum number of drug molecules to include in the ADC. Default is 0 if not specified."],
    ["Max Drug", True, "The maximum number of drug molecules to include in the ADC."],
    ["Tolerance (Da)", False, "The Tolerance in Da. Default is 50 Da if not specified."],
    ["Fixed Mod File", False, "The File Name or Path of the Fxied Mod File. "
                              "Can be either Excel or CSV. The file should have \"Mass\" and \"Name\" "
                              "columns with these exact headers. Can also have \"Number\" if a "
                              "multiple is used If not specified, no modifications will be used. "
                              "Note, it will apply one set of all fixed mods to each sequence. "
                              "If you specify, \"Seq1+Seq2\", it will apply the fixed mods to both sequences."
                              "If you specify a single amino acid sequence or float mass in the "
                              "Protein Mass column, it will apply the fixed mods only once. "],
    ["Apply Fixed Mods", False,
     "A column specifying which sequences should get the fixed mods. "
     "Should have the format of \"Seq1 Seq2 Seq3\" where Seq1 is the Sequence 1 "
     "column name, etc. Delimiters do not matter. "
     "You can specify \"All\" to apply mods to all sequences "
     "or \"None\" to apply to none. "
     "It will assume yes to all if this column is not present. "],
    ["Global Fixed Mod", False, "A column specifying a global fixed mass shift to apply to all complexes. "
                                "Unlike Fixed Mod File, which applies a fixed modification to each sequence, "
                                "this is applied only once to each complex. "
                                "Also, it is a single float value rather than a file. "
                                "This will be have identically to the Fixed Mod File "
                                "if you specify a single amino acid sequence or mass float in the "
                                "Protein Mass column. "],
    ["Sequence {n}", False,
     "The amino acid sequence or mass for the {n} protein or molecule. There can be multiple sequences, "
     "each as their own columns and with unique names. If it can convert to float, it will. "
     "If not, it will assume it is an amino acid sequence with 1 letter codes."],
    ["Disulfides Oxidized", False,
     "A column specifying the sequences that should be fully disulfide oxidized. "
     "Should have the format of \"Seq1 Seq2 Seq3\" where Seq1 is the Sequence 1 "
     "column name, etc. It will not work if only the amino acid sequence is given. Delimiters do not matter. "
     "It will assume no to all if this column is not present. "
     "You can specify \"All\" to oxidize all sequences "
     "or \"None\" to oxidize none. "
     "Will only work if sequences are amino acid codes with C. "
     "It will subtract one H mass for each C."],
]


def find_file(fname, folder, use_converted=True):
    # If use_converted is true, it will look for the converted file first, then the original.
    if use_converted:
        extensions = known_extensions[::-1]
    else:
        extensions = known_extensions

    if use_converted:
        # try to find a converted file first
        fnameroot = os.path.splitext(fname)[0]
        for ext in extensions:
            testpath = os.path.join(folder, fnameroot + ext)
            if os.path.exists(testpath):
                return testpath

    # Tries to find the extension of the file
    extension = os.path.splitext(fname)[1].lower()
    if extension in extensions:
        # If it finds a file extension it recoginzes, it will try to find the file with that extension
        if os.path.exists(os.path.join(folder, fname)):
            # Return the file with the folder
            return os.path.join(folder, fname)
        else:
            # Return the original file name and hope for the best
            return fname

    else:
        # If a file extension is not regonized, it will try to find the file with the extensions it knows
        for ext in extensions:
            if os.path.exists(os.path.join(folder, fname + ext)):
                return os.path.join(folder, fname + ext)
    return fname


def check_for_word_in_keys(df, word="Correct"):
    for k in df.keys():
        if word in k:
            return True
    return False


def set_param_from_row(eng, row):
    for k in row.keys():
        val = row[k]
        if isinstance(val, (float, int)):
            if math.isnan(val):
                continue
        if isinstance(val, str):
            if val == "nan":
                continue

        if "Config Peak Window" in k:
            try:
                eng.config.peakwindow = float(val)
            except Exception as e:
                print("Error setting peak window", k, val, e)
        if "Config Peak Thres" in k:
            try:
                eng.config.peakthresh = float(val)
            except Exception as e:
                print("Error setting peak threshold", k, val, e)
        if "Config Low Mass" in k:
            try:
                eng.config.masslb = float(val)
            except Exception as e:
                print("Error setting low mass", k, val, e)
        if "Config High Mass" in k:
            try:
                eng.config.massub = float(val)
            except Exception as e:
                print("Error setting high mass", k, val, e)

        if "Config Low m/z" in k or "Config Low mz" in k:
            try:
                eng.config.minmz = float(val)
            except Exception as e:
                print("Error setting low m/z", k, val, e)
        if "Config High m/z" in k or "Config High mz" in k:
            try:
                eng.config.maxmz = float(val)
            except Exception as e:
                print("Error setting high m/z", k, val, e)

        if "Config Sample Mass Every" in k:
            try:
                eng.config.massbins = float(val)
            except Exception as e:
                print("Error setting massbins", k, val, e)

        if "Config m/z Peak FWHM" in k or "Config mz Peak FWHM" in k:
            try:
                eng.config.mzsig = float(val)
            except Exception as e:
                print("Error setting massbins", k, val, e)

        if "Config m/z Peak Shape" in k or "Config mz Peak Shape" in k:
            try:
                eng.config.psfun = int(val)
            except Exception as e:
                print("Error setting massbins", k, val, e)

        # print(eng.config.maxmz, eng.config.minmz, k)
    return eng


def check_for_value(row, key, correcttypetuple):
    if key in row:
        if isinstance(row[key], correcttypetuple):
            return True
    return False


def check_for_floatable(row, key):
    try:
        val = float(row[key])
        if math.isnan(val):
            return False
        return True
    except Exception as e:
        return False


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
        self.vmodfile = None
        self.fmodfile = None
        self.vmoddf = None
        self.fmoddf = None
        self.autopw = True
        self.correct_pair_mode = False
        self.dar_mode = False
        self.time_range = None

    def run_file(self, file=None, decon=True, use_converted=True, interactive=False):
        self.top_dir = os.path.dirname(file)
        self.rundf = file_to_df(file)

        self.run_df(decon=decon, use_converted=use_converted, interactive=interactive)

    def run_df(self, df=None, decon=True, use_converted=True, interactive=False):

        # Print the data directory and start the clock
        clockstart = time.perf_counter()
        # Set the Pandas DataFrame
        if df is not None:
            self.rundf = df

        # Get into the right relative directory
        if os.path.isdir(self.top_dir):
            os.chdir(self.top_dir)

        # Set up the lists for the results
        htmlfiles = []

        # Check if the "Correct" column is in the DataFrame to start correct pair mode
        self.correct_pair_mode = check_for_word_in_keys(self.rundf, "Correct")
        # Check for the DAR mode
        self.dar_mode = check_for_word_in_keys(self.rundf, "Max Drugs")
        # if self.correct_pair_mode:
        #    self.rundf = remove_columns(self.rundf, "Height")

        # Loop through the DataFrame
        for i, row in self.rundf.iterrows():
            self.autopw = True
            self.eng.reset_config()
            path = self.get_file_path(row, use_converted=use_converted)

            # Get the time range
            self.time_range = get_time_range(row)

            # If the file exists, open it
            if os.path.isfile(path):
                print("Opening:", path)
                if not use_converted:
                    print("Refreshing")
                self.eng.open_file(path, time_range=self.time_range, refresh=not use_converted, silent=True)

                # If the config file is specified, load it
                if "Config File" in row:
                    try:
                        self.eng.load_config(row["Config File"])
                        print("Loaded Config File:", row["Config File"])
                        # If a config file is loaded, it will not use the auto peak width
                        self.autopw = False
                    except Exception as e:
                        print("Error loading config file", row["Config File"], e)

                # Set the deconvolution parameters from the DataFrame
                self.eng = set_param_from_row(self.eng, row)

                # Run the deconvolution or import the prior deconvolution results
                if decon:
                    # If the Config m/z Peak FWHM is specified, do not use the auto peak width
                    if "Config m/z Peak FWHM" in row:
                        self.autopw = not check_for_floatable(row, "Config m/z Peak FWHM")
                    print("Auto Peak Width", self.autopw)
                    self.eng.autorun(auto_peak_width=self.autopw, silent=True)
                else:
                    self.eng.unidec_imports(efficiency=False)
                    self.eng.pick_peaks()

                results_string = None

                # The First Recipe, correct pair mode
                if self.correct_pair_mode:
                    # Run correct pair mode
                    newrow = self.run_correct_pair(row)

                    # Merge the row back in the df
                    self.rundf = set_row_merge(self.rundf, newrow, [i])

                if self.dar_mode:
                    # Run DAR mode
                    newrow = self.run_dar(row)

                    # Merge the row back in the df
                    self.rundf = set_row_merge(self.rundf, newrow, [i])
                    try:
                        results_string = "The Drug-to-Antibody Ratio (DAR) is: " + str(newrow["DAR"])
                    except Exception:
                        results_string = None

                # Generate the HTML report
                outfile = self.eng.gen_html_report(open_in_browser=False, interactive=interactive,
                                                   results_string=results_string)
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
        # print(self.rundf)
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
        self.data_dir = ""
        if "Data Directory" in row:
            data_dir = row["Data Directory"]
            # print(data_dir, os.path.isdir(data_dir), os.getcwd())
            if os.path.isdir(str(data_dir)):
                self.data_dir = data_dir
        # Find the file
        outpath = find_file(file, self.data_dir, use_converted)
        # print("File Path:", outpath)
        return os.path.abspath(outpath)

    def get_tol(self, row):
        # Extract the tolerance for peak matching
        if "Tolerance (Da)" in row:
            try:
                self.tolerance = float(row["Tolerance (Da)"])
            except Exception:
                pass

    def get_mod_files(self, row):
        self.vmoddf = None
        # Get and read the mod file
        if "Variable Mod File" in row:
            self.vmodfile = row["Variable Mod File"]
            if os.path.isfile(str(self.vmodfile)):
                self.vmoddf = file_to_df(self.vmodfile)
                print("Loaded Variable Mod File: ", self.vmodfile)

        self.fmoddf = None
        # Get and read the mod file
        if "Fixed Mod File" in row:
            self.fmodfile = row["Fixed Mod File"]
            if os.path.isfile(str(self.fmodfile)):
                self.fmoddf = file_to_df(self.fmodfile)
                print("Loaded Fixed Mod File: ", self.fmodfile)

    def run_correct_pair(self, row, pks=None):
        if pks is None:
            pks = self.eng.pks
        # Get tolerance
        self.get_tol(row)
        # Get mod files
        self.get_mod_files(row)

        # Match to the correct peaks
        newrow = UPP_check_peaks(row, pks, self.tolerance, vmoddf=self.vmoddf, fmoddf=self.fmoddf)

        return newrow

    def run_dar(self, row, pks=None):
        if pks is None:
            pks = self.eng.pks

        # Extract the tolerance for peak matching
        self.get_tol(row)
        # Get the mod files
        self.get_mod_files(row)

        # Set the Fixed Mod if applicable
        fmod_mass = 0
        if self.fmoddf is not None:
            try:
                fmod_mass, fmod_label = parse_fmoddf(self.fmoddf)
                fmod_mass = float(fmod_mass)
            except Exception:
                print("Error calculating fixed mod mass.", self.fmoddf)
        # Add in the global fixed mod if applicable
        try:
            global_fixed_mod = parse_global_fixed_mod(row)
            global_fixed_mod = float(global_fixed_mod)
            fmod_mass += global_fixed_mod
        except Exception:
            print("Error calculating global fixed mod mass.", row)

        # Read the Protein Mass
        if "Protein Mass" in row:
            try:
                protein_mass = float(row["Protein Mass"])
                protein_mass += fmod_mass
            except Exception:
                try:
                    if "Seq" in row["Protein Mass"]:
                        masses, labels = calc_pairs(row, fmoddf=self.fmoddf, keywords=["Protein Mass"])
                        protein_mass = masses[0]
                    else:
                        protein_mass = calc_pep_mass(row["Protein Mass"])
                        protein_mass += fmod_mass

                except Exception:
                    print("Error in calculating the protein mass. Please specify Protein Mass in the DataFrame.")
                    return row
        else:
            print("Error: No Protein Mass Specified. Please specify Protein Mass in the DataFrame.")
            return row

        # Read the Drug Mass
        if "Drug Mass" in row:
            try:
                drug_mass = float(row["Drug Mass"])
            except Exception:
                print("Error: Drug Mass is not a number. Please specify Drug Mass in the DataFrame.")
                return row
        else:
            print("Error: No Drug Mass Specified. Please specify Drug Mass in the DataFrame.")
            return row

        # Read the min and max drugs
        if "Min Drugs" in row:
            try:
                min_drugs = int(row["Min Drugs"])
            except Exception:
                print("Error: Min Drugs is not an integer. Please specify Min Drugs in the DataFrame. Using 0.")
                min_drugs = 0
        else:
            min_drugs = 0

        if "Max Drugs" in row:
            try:
                max_drugs = int(row["Max Drugs"])
            except Exception:
                print("Error: Max Drugs is not an integer. Please specify Max Drugs in the DataFrame.")
                return row
        else:
            print("Error: No Max Drugs Specified. Please specify Max Drugs in the DataFrame.")
            return row

        if protein_mass > 0 and drug_mass > 0 and min_drugs >= 0 and max_drugs > 0:
            print("Running DAR Calculation. Min Drugs:", min_drugs, "Max Drugs:", max_drugs, "Protein Mass:",
                  protein_mass,
                  "Drug Mass:", drug_mass)
            dar_val = dar_calc(pks, protein_mass, drug_mass, min_drugs, max_drugs, self.tolerance)
            row["DAR"] = dar_val
            return row
        else:
            print("Error: DAR Calculation Failed. Please check your parameters.")
            print("Min Drugs:", min_drugs, "Max Drugs:", max_drugs, "Protein Mass:", protein_mass,
                  "Drug Mass:", drug_mass)
            return row


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
        path = "C:\\Data\\Wilson_Genentech\\sequences_short3.xlsx"
        path = "C:\\Data\\Wilson_Genentech\\BsAb\\BsAb test short.xlsx"
        # path = "C:\\Data\\Wilson_Genentech\\BsAb\\BsAb test.xlsx"
        # path = "C:\\Data\\Wilson_Genentech\\BsAb\\test2.csv"
        path = "C:\\Data\\Wilson_Genentech\\Test\\BsAb test v2.xlsx"
        path = "C:\\Data\\UPPDemo\\BsAb\\BsAb test short.xlsx"
        path = "C:\\Data\\Wilson_Genentech\\DAR\\Biotin UPP template test.xlsx"
        pd.set_option('display.max_columns', None)
        batch.run_file(path, decon=True, use_converted=True, interactive=False)

        batch.open_all_html()
        pass
