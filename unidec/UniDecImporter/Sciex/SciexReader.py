import clr
import os
import sys

curr_dir = os.path.realpath(__file__)

if hasattr(sys, 'frozen') and hasattr(sys, '_MEIPASS'):
    curr_dir = curr_dir.split("\\")
    curr_dir = "\\".join(curr_dir[:-4])
    # print("Running in PyInstaller at", curr_dir)
else:
    curr_dir = curr_dir.split("\\")
    curr_dir = "\\".join(curr_dir[:-1])
    # print("Running in Python at", curr_dir)


dll_files = [
    "Clearcore2.Muni",
    "BlaisWiff",
    "Clearcore2.Data.AnalystDataProvider",
    "Clearcore2.Data",
    "Clearcore2.Utility",
    "Clearcore2.RawXYProcessing",
    "Clearcore2.StructuredStorage",
    "Sciex.TofTof.T2DFMan",
    "WiffReaderCOM",
]

try:
    for dll in dll_files:
        dll_path = os.path.join(curr_dir, dll) + ".dll"
        clr.AddReference(dll_path)
except Exception as e:
    print("Failed to Sciex load DLLs")
    print(e)


import System

try:
    from WiffReader import BlaisModel
except Exception as e:
    print("Failed to import BlaisModel")
    print(e)


class SciexReader():
    def __init__(self, path, selected_sample=None, selected_experiment=None):
        print("Reading Sciex Data:", path)
        self.reader = None
        try:
            self.reader = BlaisModel()
            self.reader.OpenWiffFile(path)
        except Exception as e:
            print(f"Error while opening WIFF file: {e}")
        self.num_samples = self.reader.GetSamples()
        self.sampleTypes = self.populate_sample_dict()
        self.scan_dict = self.get_curr_scan_dict(self.sampleTypes)
        self.all_scanranges = self.get_all_scan_dict()

    def populate_sample_dict(self):
        sample_dict = {}

        total_samples = self.reader.GetSamples()

        for sample in range(total_samples):
            num_exp = self.reader.GetExperiments(sample)
            sample_dict[sample] = list(range(num_exp))  # Store all experiment indices for the sample

        return sample_dict

    def get_curr_scan_dict(self, example_dict):
        scan_dict = {}
        for ex in example_dict.keys():
            scan_dict[ex] = self.reader.GetNumCycles(ex)
        return scan_dict

    def get_all_scan_dict(self):
        """Stores { (Sample Number, Experiment Number) : Scan Range } in a dictionary."""
        scan_info_dict = {}

        #print("Total Samples:", self.num_samples)

        for sample, experiments in self.sampleTypes.items():
            #use this line to see ALL samples and experiments from the .wiif file
            #print(f"Sample(s) : {sample + 1}:  Experiment(s): {len(experiments)}")

            for experiment in experiments:
                scan_range = [0, self.reader.GetNumCycles(sample)]
                scan_info_dict[(sample + 1, experiment + 1)] = {
                    "Scan Range": scan_range
                }
        return scan_info_dict


if __name__ == '__main__':
    path = "C:\\Users\\MartyLabsOfficePC\\OneDrive - University of Arizona\\Desktop\\20230816_Myoglobin 0666 ugmL 01.wiff"
    sciex = SciexReader(path)
    sampletypes = sciex.reader.GetSampleNames()
    exp_info = sciex.reader.ExperimentInfo(0, 0)
    for i in exp_info:
        print(i)

