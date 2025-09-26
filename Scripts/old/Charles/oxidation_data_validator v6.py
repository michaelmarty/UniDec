import pandas as pd
import numpy as np
import pprint
from matplotlib import pyplot as plt
import wx
import os
import datetime
from pyteomics import fasta
from Bio.PDB import *
from copy import deepcopy

''' v2 has a class PercentOxidationAnalyzer that can be called from other files and two seperate fucntions:
        1) OxidationValidatorSingle    -  Analyzes for specific protein and returns list of percent oxidation and the oxidized amino acid
        2) OxidationValidator          -  Identifies Unique proteins and analyzes for them all

    v3 has new class OxidizedAminoAcidGrapher(wx.Frame) with one method:
        1) BarGrapher   (Doesnt Work)  -  Creates a bar graph with oxidation on y and oxidized amino acid
    
       ... as well as a new function:
        1) SaveToExcel                 -  Saves a chart of [Oxidized Amino Acid : Percent Oxidation] to excelsheet. 
                                          If no specified file, it saves to desktop
                                          
    v4 attempts to solve the in-acuracy that comes from not acounting for the same amino acids from different overlapping segments 
    
    v5 fixes the SaveToExcel method and corrects algorithm calculating percecnt modified and makes Oxidation Validator single return a list of oxidation values paired with segments
    
    v6 fixes BarGrapher from v3 and renames two methods.
       
       BarGrapher() FIXED
      
       Two renamed methods:
       1) OxidationValidatorSingle()  --> AminoAcidOxidationSingle()
       2) OxidationValidator()        --> AminoAcidOxitation
    
       MTM
       
       Added substring method for matching protein names
       Changed bar plots to read in sequence and plot it
       Added fasta import for global protien sequence data
       Added PDB import, write %ox to beta factor, and export
'''

protection_factors={"C":29.2, "M":20.5, "W":17.4, "Y":12.0, "F":11.2, "H":9.3, "L":4.4, "I":4.4, "R":2.9, "K":2.2
                    , "V":1.9, "T":1.6, "S":1.4, "P":1.0, "E":0.69, "Q":0.66, "N":0.44, "D":0.42, "A":0.14, "G":0.04}


def substring(array, string):
    array = np.array(array,dtype=str)
    return np.flatnonzero(np.core.defchararray.find(array, string) != -1)

def safedivide(a, b):
    """
    Divide numpy arrays b from a. Define a/0=0 to avoid zero division error.
    :param a: Numerator
    :param b: Denominator
    :return: a/b
    """
    c = deepcopy(b)
    c[b != 0] = a[b != 0] / b[b != 0]
    return c

class PercentOxidationAnalyzer:
    def __init__(self, filePath, exclusionTimeDif=0.25):
        self.filePath = filePath
        self.exclusionTimeDif = exclusionTimeDif
        print("File accessed from " + self.filePath + "\n")

    def OxidationValidatorSingle(self, protein_name, exact = False):
        EXCLUSION_TIME_DIFFERENCE = self.exclusionTimeDif
        self.protein_name = protein_name

        # Lab file location.
        # excel_file = 'Z:Group Share\\James Rohrbough\\FPOP\\MSMS raw data\\ubiquitin\\FPOP peptides\\JGR_2319383_2Ox-calib_Peptides.xlsx'

        peptides = pd.read_excel(self.filePath)

        # sorted_by_protein = peptides.sort(['Protein Name'])
        # sorted_by_RT = peptides.sort.value(['Scan Retention'])

        # Grab filtered data.
        if exact:
            filtered_data = peptides.loc[peptides["Protein Name"] == self.protein_name]
        else:
            filtered_data = peptides.loc[substring(peptides["Protein Name"], self.protein_name)]
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

        # Initialize lists
        nresidues = 76 #TODO: find this value automatically
        oxidized_amino_acids_ion_currents = np.ndarray(shape=nresidues)
        oxidized_amino_acids_ion_currents.fill(0)
        total_sum = np.ndarray(shape=nresidues)
        total_sum.fill(0)

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
            print("'{}' reference retention time: {}".format(key, reference_retention_times[key]))

            for index, row in group.iterrows():
                # Grab start residue for current peptide.
                frag_len = row["Start and End Residues In Protein"].replace("[", "").replace("]", "").split(" ")
                start = int(frag_len[0])-1
                end = int(frag_len[2])-1

                if abs(row["Scan Retention Time"] - reference_retention_times[key]) > EXCLUSION_TIME_DIFFERENCE:
                    essential_sequence = row["Essential Sequence"]
                    parsed_sequence = essential_sequence

                    # Include modified fragments only.
                    if pd.isnull(row["Mods"]) is False:
                        frag_len = start
                        total_ion_current = row["Total Ion Current"]
                        # Find oxidation based on brackets in column "Total Ion Current".
                        while "]" in parsed_sequence:
                            parsed_sequence = parsed_sequence.split("[", 1)
                            frag_len += len(parsed_sequence[0])
                            parsed_sequence = parsed_sequence[1].split("]", 1)
                            parsed_sequence = parsed_sequence[1]
                            oxidized_amino_acids_ion_currents[frag_len-1] += total_ion_current
                            #print(essential_sequence, frag_len, total_ion_current)

                        total_sum[start:end] += total_ion_current #TODO: Note rewrite here
                    else:
                        # TODO: Check MTM Rewrite
                        total_sum[start:end] += row["Total Ion Current"]

                else:
                    #TODO: Check MTM Rewrite
                    #print(row["Essential Sequence"], start, end, row["Total Ion Current"])
                    total_sum[start:end] += row["Total Ion Current"]

                #if start<67 and end>67:
                #    print(row["Essential Sequence"], start, end, row["Total Ion Current"], total_sum[38])


        #print(total_sum[38], oxidized_amino_acids_ion_currents[38])
        percent_oxidation = safedivide(oxidized_amino_acids_ion_currents,total_sum) # Added safe divide for parts that are never hit
        #TODO: how do we handle sequence parts that aren't covered?
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




def BarGrapher(protein_name, oxidized_list, seq= None):
    if seq is None:
        plt.bar(range(oxidized_list.size), oxidized_list)
    else:
        slist = [s + "\n" +str(i+1) for i, s in enumerate(seq)]
        plt.bar(range(oxidized_list.size), oxidized_list, tick_label = slist )

    plt.xlabel('Amino Acid')
    plt.ylabel('Percent Modified')
    plt.title('Percent Oxidized by Amino Acid in ' + protein_name)
    plt.suptitle(protein_name)

    plt.show()



def SaveToExcel(protein_name, oxidized_amino_acid_dict):

    oxidized_amino_acid_dict #= oxidized_amino_acid_dict.items()
    df = pd.DataFrame(oxidized_amino_acid_dict)
    df = df.rename(columns={0:'Percent Oxidized', 1:'Amino Acid Location'})
    df = df.sort_values(by=['Percent Oxidized'])

    date = datetime.datetime.now()
    writer = pd.ExcelWriter(path="/Users/charles/Desktop/Polyubiquitin-B|.xlsx", date_format=date, engine='xlsxwriter', datetime_format='YYYY-MM-DD', mode='w')
    #protein_name = protein_name + ".xlsx"
    df.to_excel(protein_name+".xlsx")
    df.to_excel(writer, sheet_name=protein_name)
    print(df)

def import_fasta(file_name):
    fd = fasta.read(file_name)
    return(fd)

def get_sequence_from_fasta(file_name, protein_name):
    fd = import_fasta(file_name)
    for f in fd:
        des = f.description
        if protein_name in des:
            seq = f.sequence
            print(seq)
    return seq

def protein_factors_from_seq(sequence):
    pf = []
    for l in sequence:
        pf.append(protection_factors[l])
    return np.array(pf)

def write_to_pdb(pdbfile, array,  name="UBI"):
    parser = PDBParser()
    structure = parser.get_structure(name, pdbfile)
    factor = 100./np.amax(array)
    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain):
                if i<len(array):
                    for j, atom in enumerate(residue):
                        atom.set_bfactor(array[i]*factor)

    newname= os.path.splitext(pdbfile)[0]+"_fpop.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(newname)

def main():
    #app = wx.App()

    #app.MainLoop()



    '''Analysis'''
    protein_name = "Polyubiquitin-B"

    seq = get_sequence_from_fasta("ubiquitin.FASTA", protein_name)
    pf = protein_factors_from_seq(seq)
    Analyzer = PercentOxidationAnalyzer("output.xlsx")
    #Analyzer.OxidationValidator()
    #oxidized_aminoacids = dict()

    oxidized_aminoacids = Analyzer.OxidationValidatorSingle(protein_name)
    #print(pprint.pformat(oxidized_aminoacids))
    oxidized_aminoacids = oxidized_aminoacids / pf

    write_to_pdb("1ubq.pdb", oxidized_aminoacids)

    '''Grapher'''
    BarGrapher(protein_name, oxidized_aminoacids, seq)

    '''Save To Excel'''
    #SaveToExcel(protein_name, oxidized_aminoacids)
    #print(len(oxidized_aminoacids))

if __name__ == "__main__":
    main()
