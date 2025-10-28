import pandas as pd
import numpy as np
import pprint
import pyautogui
import wx
import os
import datetime
import xlsxwriter

''' v2 has a class PercentOxidationAnalyzer that can be called from other files and two seperate fucntions:
        1) OxidationValidatorSingle -  analyzes for specific protein and returns list of percent oxidation and the oxidized amino acid
        2) OxidationValidator       -  Identifies Unique proteins and analyzes for them all

    v3 has new class OxidizedAminoAcidGrapher(wx.Frame) wiht one method:
        1) BarGrapher   (Doesnt Work)               -  Creates a bar graph with oxidation on y and oxidized amino acid
    
    ... as well as a new function:
        1) SaveToExcel                 -  Saves a chart of [Oxidized Amino Acid : Percent Oxidation] to excelsheet. 
                                          If no specified file, it saves to desktop
    v4 attempts to solve the in-acuracy that comes from not acounting for the same amino acids from different overlapping segments 
'''


class PercentOxidationAnalyzer:
    def __init__(self, filePath, exclusionTimeDif):
        self.filePath = filePath
        self.exclusionTimeDif = exclusionTimeDif
        print("File accessed from " + self.filePath + "\n")

    def OxidationValidatorSingle(self, protein_name):
        '''EXCLUSION_TIME_DIFFERENCE = 0.25'''
        EXCLUSION_TIME_DIFFERENCE = self.exclusionTimeDif
        self.protein_name = protein_name

        # Lab file location.
        # excel_file = 'Z:Group Share\\James Rohrbough\\FPOP\\MSMS raw data\\ubiquitin\\FPOP peptides\\JGR_2319383_2Ox-calib_Peptides.xlsx'
        # Local file location.

        '''filePath = "/Users/charles/PycharmProjects/MSProteinOxidation/output.xlsx"'''
        peptides = pd.read_excel(self.filePath)

        # sorted_by_protein = peptides.sort(['Protein Name'])
        # sorted_by_RT = peptides.sort.value(['Scan Retention'])

        # Grab filtered data.
        filtered_data = peptides.loc[peptides["Protein Name"] == self.protein_name]
        '''filtered_data = peptides.loc[peptides["Protein Name"] == "Polyubiquitin-B|"]'''
        filtered_data = filtered_data.fillna("")

        # Summarize results.
        print("{} peptides extracted from a total of {}".format(
            len(filtered_data), len(peptides)
        ) +"\n")

        # Store the counts for each base sequence. Maps sequence to count.
        base_sequence_counts = dict()

        for base_sequence in filtered_data["Base Sequence"]:
            # Create dictionary entry if it doesn't exist.
            base_sequence_counts.setdefault(base_sequence, 0)
            # Update the count.
            base_sequence_counts[base_sequence] += 1

        # Maps the base_sequence to the retention time.
        reference_retention_times = dict()

        # Initialize lists that will be used to calculate percent modification
        oxidized_amino_acids_ion_currents = np.ndarray(shape=76+1)      #76 is the amount of amino acids in Polyubiquitin-B
        oxidized_amino_acids_ion_currents.fill(0)
        total_sum = 0.0



        for key in base_sequence_counts:
            # Get group of matching base sequences
            group = peptides.loc[peptides["Base Sequence"] == key]


            # Find the retention times for each base sequence.
            for index, row in group.iterrows():
                # Append unmodified entries to calculate reference retention time.
                if pd.isnull(row["Mods"]):
                    reference_retention_times.setdefault(key, list())
                    reference_retention_times[key].append(row["Scan Retention Time"])

            # Error Handling: When a reference time couldn't be found.
            if reference_retention_times.get(key) is None:
                print(
                    "EXCLUSION: Could not find base sequence '{}' without"
                    " mods to use as reference time.".format(key))
                continue

            # Reference retention time = mean(unmodified_entries)
            reference_retention_times[key] = np.mean(reference_retention_times[key])

            # Notify user of retention time.
            print("'{}' reference retention time: {}".format(
                key, reference_retention_times[key]))



            for index, row in group.iterrows():
                if (
                    row["Scan Retention Time"] <
                    (reference_retention_times[key] - EXCLUSION_TIME_DIFFERENCE) or
                    row["Scan Retention Time"] >
                    (reference_retention_times[key] + EXCLUSION_TIME_DIFFERENCE)
                ):
                    essential_sequence = row["Essential Sequence"]
                    parsed_sequence = essential_sequence

                    # Grab start residue for current peptide.
                    frag_len = int(
                        row["Start and End Residues In Protein"].split(" ")[0][1:]) - 1

                    # Include modified fragments only.
                    if pd.isnull(row["Mods"]) is False:
                        # Find oxidation based on brackets in column "Total Ion Current".
                        while "]" in parsed_sequence:
                            total_ion_current = row["Total Ion Current"]
                            print("total ion current: " + str(total_ion_current))
                            total_sum += total_ion_current
                            print("total sum: " + str(total_sum))

                            parsed_sequence = parsed_sequence.split("[", 1)
                            frag_len += len(parsed_sequence[0])
                            parsed_sequence = parsed_sequence[1].split("]", 1)
                            parsed_sequence = parsed_sequence[1]

                            oxidized_amino_acids_ion_currents[frag_len] += total_ion_current


        percent_oxidation = oxidized_amino_acids_ion_currents/total_sum
        #value = np.sum(percent_oxidation)
        print(percent_oxidation)
        return percent_oxidation

    def OxidationValidator(self):
        peptides = pd.read_excel(self.filePath)
        protein_list = peptides["Protein Name"]
        #print(protein_list)
        protein_list = protein_list.astype(str)
        protein_list = np.unique(protein_list)

        print("List of unique proteins: " + protein_list)

        for protein in protein_list:
            if protein == "nan":
                pass
            else:
                print("\n" "Currently Analyzing Protein " + protein + "\n")
                self.oxidationValidatorSingle()

#class OxidizedAminoAcidGrapher(wx.Frame):

    #def __init__(self):
    #def Grapher(wx.Frame)

def SaveToExcel(protein_name, oxidized_amino_acid_dict):

    oxidized_amino_acid_dict = oxidized_amino_acid_dict.items()
    df = pd.DataFrame(oxidized_amino_acid_dict)
    df = df.rename(columns={0:'Amino Acid Location', 1:"Percent Oxidized"})
    df = df.sort_values(by=['Amino Acid Location'])

    date = datetime.datetime.now()
    writer = pd.ExcelWriter(path="/Users/charles/Desktop/", date_format=date, engine='xlsxwriter', datetime_format='YYYY-MM-DD', mode='w')

    df.to_excel(excel_writer=writer, sheet_name=protein_name)
    #print(df)

def main():

    app = wx.App()
    #grapher = Grapher(None)
    #grapher.Show()
    app.MainLoop()

    protein_name = "Polyubiquitin-B|"
    Analyzer = PercentOxidationAnalyzer("/Users/charles/PycharmProjects/MSProteinOxidation/output.xlsx", 0.25)

    #Analyzer.OxidationValidator()

    oxidized_aminoacids = dict()
    oxidized_aminoacids = Analyzer.OxidationValidatorSingle(protein_name)

    #print(pprint.pformat(oxidized_aminoacids))
    #SaveToExcel(protein_name, oxidized_aminoacids)
    #print(len(oxidized_aminoacids))

if __name__ == "__main__":
    main()
