import pandas as pd
import numpy as np
import pprint
from matplotlib import pyplot as plt
from operator import itemgetter
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
    
    v7
        MTM
       
       Added substring method for matching protein names
       Changed bar plots to read in sequence and plot it
       Added fasta import for global protien sequence data
       Added PDB import, write %ox to beta factor, and export
       
    v8  
       
       BarGrapher() can now take dictionaries and graph them
       Added Peptide analyzer (v2 of that method)
       
    v9 
       Adapts PeptideValidator() to be able to read excel files with different header names 
       
    v10
       Simplifies/Improves readability of PeptideValidator() by creating external methods
       Combines analysis of peptide and analysis of amino acid to both execute when PeptideValidator() is called

    v11 
        Deletes redundant AminoOxidationValidator() method from v10
        Continue to simplify PecentOxidationAnalyzer's code
'''

protection_factors={"C":29.2, "M":20.5, "W":17.4, "Y":12.0, "F":11.2, "H":9.3, "L":4.4, "I":4.4, "R":2.9, "K":2.2
                    , "V":1.9, "T":1.6, "S":1.4, "P":1.0, "E":0.69, "Q":0.66, "N":0.44, "D":0.42, "A":0.14, "G":0.04}


def substring(array, string):
    array = np.array(array, dtype=str)
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
    def __init__(self, filePath, exclusionTimeDif, FASTA):
        self.filePath = filePath
        self.exclusionTimeDif = exclusionTimeDif
        print("File accessed from " + self.filePath + "\n")
        self.FASTA = FASTA

    def PeptideValidator(self, protein_name, exact=False):
        self.exact = exact
        self.protein_name = protein_name
        areDUPLICATES = False
        analyzeAMINO = True
        # Grab full seq from imported Fasta file
        self.sequence = get_sequence_from_fasta(self.FASTA, protein_name)
        seq = self.sequence
        # Lab file location.
        # excel_file = 'Z:Group Share\\James Rohrbough\\FPOP\\MSMS raw data\\ubiquitin\\FPOP peptides\\JGR_2319383_2Ox-calib_Peptides.xlsx'
        peptides = pd.read_excel(self.filePath)
        # sorted_by_protein = peptides.sort(['Protein Name'])
        # sorted_by_RT = peptides.sort.value(['Scan Retention'])

        # Identify and set column header names
        self.get_column_header_names(peptides)
        header_base_sequence = self.header_base_sequence
        header_retention_time = self.header_retention_time
        header_modifications  = self.header_modifications

        #protein_accession = "P0CG53"       #Polyubiquitin-Bov
        protein_accession = "P62894"        #CytC-Bov      '''Either find a way to get master protein accession from fasta or demand it as an input along with/instead of protein name'''
        desired_proteins = [self.protein_name, protein_accession]

        # Cuts dataset to only include data from desired protein
        filtered_data = self.select_proteins(desired_proteins, peptides)

        # Summarizes how many peptides are being extracted from data set
        print("{} peptides extracted from a total of {}".format(
            len(filtered_data), len(peptides)) +"\n")


        # Maps the base_sequence to the retention time.
        reference_retention_times = dict()

        # Assigns key to each unique peptide sequence
        unique_peptide_list = dict()
        peptide_list = filtered_data[header_base_sequence]
        peptide_list = np.unique(peptide_list)
        # Sort by starting position
        self.peptide_list = self.find_loc_and_sort(peptide_list, seq)
        print(peptide_list)
        key = 0
        for base_sequence in self.peptide_list:

            unique_peptide_list.setdefault(base_sequence, key)
            key += 1

        #Count how many times a peptide occors
        base_sequence_counts = dict()
        for base_sequence in filtered_data[header_base_sequence]:
            for spot, pep_name in self.peptide_list.items():
                if(pep_name) == base_sequence:
                    # Create dictionary entry if it doesn't exist.
                    base_sequence_counts.setdefault(spot, 0)
                    # Update the count.
                    base_sequence_counts[spot] += 1

        # Initialize lists for peptide
        npeptides = len(unique_peptide_list)
        self.oxidized_peptide_ion_currents = np.ndarray(shape=npeptides)
        self.oxidized_peptide_ion_currents.fill(0)
        self.total_sum = 0

        # Initialize lists for amino analysis
        nresidues = len(seq)
        self.oxidized_amino_acids_ion_currents = np.ndarray(shape=nresidues)
        self.oxidized_amino_acids_ion_currents.fill(0)
        self.total_amino_sum = np.ndarray(shape=nresidues)
        self.total_amino_sum.fill(0)


        for key in base_sequence_counts:
            # Get group of matching base sequences
            group = self.eliminate_duplicates(peptides, key, ELIM_DUP = areDUPLICATES)
            print(key)
            # Find the retention times for each base sequence.
            for index, row in group.iterrows():
                # Append unmodified entries to calculate reference retention time.
                if pd.isnull(row[header_modifications]):
                    reference_retention_times.setdefault(key, list())
                    reference_retention_times[key].append(row[header_retention_time])

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

            # Create a variable that represents the current peptides spot in the dictionary
            peptide_list_key = unique_peptide_list[key]

            # Extract ion currents from excel sheet and assign them to the peptide
            self.assign_pep_ion_currents(specific_retention_times=reference_retention_times[key], peptide_list_key=peptide_list_key, peptide_group=group)

            # If selected for, also analyze percent modified by amino acid
            if(analyzeAMINO):
                self.assign_amino_ion_currents(specific_retention_times=reference_retention_times[key], peptide_group=group)

        # Calculate percent modified for each peptide
        percent_oxidation = self.oxidized_peptide_ion_currents/self.total_sum

        # Pair peptide with its oxidation value
        percent_oxidation_dict = dict()
        counter = 0
        for peptide in unique_peptide_list:
            percent_oxidation_dict.setdefault(peptide, 0)
            percent_oxidation_dict[peptide] = percent_oxidation[counter]
            counter += 1

        # Notify user of peptide oxidation values
        print(percent_oxidation_dict)

        # Graph percent modified by peptide
        BarGrapher(self.protein_name, percent_oxidation_dict)

        ##################################################

        # Calculate percent modified for each amino acid
        percent_oxidation_amino = safedivide(self.oxidized_amino_acids_ion_currents,self.total_amino_sum) # Added safe divide for parts that are never hit
        #pf = protein_factors_from_seq(seq)
        #percent_oxidation_amino = percent_oxidation_amino/pf

        # Pair amino acid with its oxidation value

        percent_oxidation_amino_dict = dict()
        slist = [s + "\n" +str(i+1) for i, s in enumerate(seq)]
        counter = 0
        for s in slist:
            percent_oxidation_amino_dict.setdefault(s, 0)
            percent_oxidation_amino_dict[s] = percent_oxidation_amino[counter]
            counter += 1

        #Notify user of amino acid oxidation values

        print("Dictionarry " + str(percent_oxidation_amino_dict))

        # Graph percent modified by amino acid
        BarGrapher(self.protein_name, percent_oxidation_amino_dict)

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

    def check_retention_time(self, specific_reference_retention_times, header_retention_time, peptide):
        if abs(peptide[header_retention_time] - specific_reference_retention_times) > self.exclusionTimeDif:
            return True
        else:
            return False

    def get_column_header_names(self, peptides):
        if 'Protein Name' in peptides.columns:
            self.header_protein_name = "Protein Name"
        elif 'Master Protein Accessions' in peptides.columns:
            self.header_protein_name = "Master Protein Accessions"
        if 'Base Sequence' in peptides.columns:
            self.header_base_sequence = "Base Sequence"
        elif 'Sequence' in peptides.columns:
            self.header_base_sequence = "Sequence"
        if 'Scan Retention Time' in peptides.columns:
            self.header_retention_time = "Scan Retention Time"
        elif 'RT [min]' in peptides.columns:
            self.header_retention_time = "RT [min]"
        if 'Mods' in peptides.columns:
            self.header_modifications = "Mods"
        elif 'Modifications' in peptides.columns:
            self.header_modifications = "Modifications"
        if 'Total Ion Current' in peptides.columns:
            self.header_ion_current = "Total Ion Current"
        elif 'Precursor Abundance' in peptides.columns:
            self.header_ion_current = "Precursor Abundance"
        if 'Full Sequence' in peptides.columns:
            self.header_full_sequence = "Full Sequence"
        else:
            self.header_full_sequence = "NULL"
            print("File has no header named: Full Sequence")

    def assign_pep_ion_currents(self, specific_retention_times, peptide_list_key, peptide_group):
        if self.header_full_sequence == "Full Sequence":
                for index, row in peptide_group.iterrows():
                    # Grab start residue for current peptides
                    if self.check_retention_time(specific_reference_retention_times=specific_retention_times, header_retention_time=self.header_retention_time, peptide=row):
                        essential_sequence = row[self.header_full_sequence]
                        parsed_sequence = essential_sequence

                        # Include modified fragments only.
                        if pd.isnull(row[self.header_modifications]) is False:
                            total_ion_current = row[self.header_ion_current]
                            # Find oxidation based on brackets in column "Total Ion Current".
                            while "]" in parsed_sequence:
                                parsed_sequence = parsed_sequence.split("[", 1)
                                parsed_sequence = parsed_sequence[1].split("]", 1)
                                parsed_sequence = parsed_sequence[1]
                                self.oxidized_peptide_ion_currents[peptide_list_key] += total_ion_current
                                #print(essential_sequence, frag_len, total_ion_current)

                            self.total_sum += total_ion_current #TODO: Note rewrite here
                        else:
                            # TODO: Check MTM Rewrite
                            self.total_sum += row[self.header_ion_current]

                    else:
                        #TODO: Check MTM Rewrite
                        #print(row["Essential Sequence"], start, end, row["Total Ion Current"])

                        self.total_sum += row[self.header_ion_current]
        elif self.header_full_sequence == "NULL":
                for index, row in peptide_group.iterrows():
                    # Grab start residue for current peptides
                    if self.check_retention_time(specific_reference_retention_times=specific_retention_times, header_retention_time=self.header_retention_time, peptide=row):
                        essential_modifications = row[self.header_modifications]

                        # Include modified fragments only.
                        if pd.isnull(row[self.header_modifications]) is False:
                            total_ion_current = row[self.header_ion_current]
                            # Find oxidation based on number of listed oxidations in column "Modifications".
                            oxidation_list = essential_modifications.split("; ")
                            for oxidation in oxidation_list:
                                oxidation_type = oxidation.split("(")
                                if "Tri" in oxidation_type[1]:
                                    self.oxidized_peptide_ion_currents[peptide_list_key] += total_ion_current*3
                                elif "Di" in oxidation_type[1]:
                                    self.oxidized_peptide_ion_currents[peptide_list_key] += total_ion_current*2
                                else:
                                    self.oxidized_peptide_ion_currents[peptide_list_key] += total_ion_current

                            self.total_sum += total_ion_current #TODO: Note rewrite here
                        else:
                            # TODO: Check MTM Rewrite
                            self.total_sum += row[self.header_ion_current]

                    else:
                        #TODO: Check MTM Rewrite
                        #print(row["Essential Sequence"], start, end, row["Total Ion Current"])
                        print("RETENTION TIME NOT MET")
                        self.total_sum += row[self.header_ion_current]

    def assign_amino_ion_currents(self, specific_retention_times, peptide_group):
        if self.header_full_sequence == "Full Sequence":
            for index, row in peptide_group.iterrows():
                # Grab start residue for current peptide.
                frag_len = row["Start and End Residues In Protein"].replace("[", "").replace("]", "").split(" ")
                print(row[self.header_protein_name])
                start = int(frag_len[0])-1
                end = int(frag_len[2])-1

                if self.check_retention_time(specific_reference_retention_times=specific_retention_times, header_retention_time=self.header_retention_time, peptide=row):
                    essential_sequence = row[self.header_full_sequence]
                    parsed_sequence = essential_sequence

                    # Include modified fragments only.
                    if pd.isnull(row[self.header_modifications]) is False:
                        frag_len = start
                        total_ion_current = row[self.header_ion_current]
                        # Find oxidation based on brackets in column "Total Ion Current".
                        while "]" in parsed_sequence:
                            parsed_sequence = parsed_sequence.split("[", 1)
                            frag_len += len(parsed_sequence[0])
                            parsed_sequence = parsed_sequence[1].split("]", 1)
                            parsed_sequence = parsed_sequence[1]
                            self.oxidized_amino_acids_ion_currents[frag_len-1] += total_ion_current
                            #print(essential_sequence, frag_len, total_ion_current)

                        self.total_amino_sum[start:end] += total_ion_current #TODO: Note rewrite here
                    else:
                        # TODO: Check MTM Rewrite
                        self.total_amino_sum[start:end] += row[self.header_ion_current]

                else:
                    #TODO: Check MTM Rewrite
                    #print(row["Essential Sequence"], start, end, row["Total Ion Current"])
                    self.total_amino_sum[start:end] += row[self.header_ion_current]

        elif self.header_full_sequence == "NULL":
            for index, row in peptide_group.iterrows():
                # Grab start residue for current peptide.
                frag_len = self.find_loc_single(peptide=row[self.header_base_sequence], sequence=self.sequence)

                if self.check_retention_time(specific_reference_retention_times=specific_retention_times, header_retention_time=self.header_retention_time, peptide=row):
                    essential_modifications = row[self.header_modifications]
                    # Include modified fragments only.
                    if pd.isnull(row[self.header_modifications]) is False:
                        total_ion_current = row[self.header_ion_current]
                        # Find oxidation based on number of listed oxidations in column "Modifications".
                        oxidation_list = essential_modifications.split("; ")
                        if len(oxidation_list) == 1:
                            for oxidation in oxidation_list:
                                start = int(frag_len[0])
                                end = int(frag_len[1])
                                oxidation_type = oxidation.split("(")
                                if "Tri" in oxidation_type[1]:
                                    pass
                                elif "Di" in oxidation_type[1]:
                                    pass
                                else:
                                    while start <= end:
                                        self.oxidized_amino_acids_ion_currents[start] += total_ion_current
                                        start += 1
                        else:
                            pass
                        while start <= end:
                            self.total_amino_sum[start] += total_ion_current #TODO: Note rewrite here
                            start += 1
                    else:
                        # TODO: Check MTM Rewrite
                        start = int(frag_len[0])
                        end = int(frag_len[1])
                        while start <= end:
                            self.total_amino_sum[start] += row[self.header_ion_current] #TODO: Note rewrite here
                            start += 1

                else:
                    #TODO: Check MTM Rewrite
                    #print(row["Essential Sequence"], start, end, row["Total Ion Current"])
                    print("RETENTION TIME NOT MET")
                    start = int(frag_len[0])
                    end = int(frag_len[1])
                    while start <= end:
                            self.total_amino_sum[start] += row[self.header_ion_current] #TODO: Note rewrite here
                            start += 1

    def eliminate_duplicates(self, peptides, key, ELIM_DUP):

        group = peptides.loc[peptides[self.header_base_sequence] == self.peptide_list[key]]
        if(ELIM_DUP):
            group = group.drop_duplicates(subset=[self.header_retention_time, self.header_ion_current], keep='first')
        return group

    def select_proteins(self, desired_protein, peptides):
        # desired_protien is a list() where the inner list contains two index. Index 0 => Protein name, Index 1 => Master Protein Accession
        protein_index = 0
        accession_index = 1
        if 'Protein Name' in self.header_protein_name:
            #Extract Protein names
            if self.exact:
                filtered_data = peptides.loc[peptides[self.header_protein_name] == desired_protein[protein_index]]
            else:
                filtered_data = peptides.loc[substring(peptides[self.header_protein_name], desired_protein[protein_index])]
                filtered_data = filtered_data.fillna("")
        elif 'Accessions' in self.header_protein_name:
            if self.exact:
                filtered_data = peptides.loc[peptides[self.header_protein_name] == desired_protein[accession_index]]
            else:
                filtered_data = peptides.loc[substring(peptides[self.header_protein_name], desired_protein[accession_index])]
                filtered_data = filtered_data.fillna("")
        else:
            print("ERROR: No column label that can be used to identify unique protein")
        return filtered_data

    def find_loc_single(self, peptide, sequence):
        starting_loc = sequence.find(peptide) + 1
        ending_loc = starting_loc + len(peptide) - 1
        peptide_location = [starting_loc, ending_loc]

        return peptide_location

    def find_loc_and_sort(self, unique_peptide_list, seq):
        LABELED_unique_peptide_list = list()
        SORTED_unique_peptide_list = list()
        SORTED_peptide_dict = dict()
        key = 0

        for peptide in unique_peptide_list:
            starting_loc = seq.find(peptide) + 1
            ending_loc = starting_loc + len(peptide) - 1
            peptide_location = [starting_loc, ending_loc]

            LABELED_unique_peptide_list.append(peptide_location)
            key += 1

        SORTED_unique_peptide_list = sorted(LABELED_unique_peptide_list, key=itemgetter(0,1))
        counter = 0
        for peptide in SORTED_unique_peptide_list:
            starting_loc = peptide[0]
            ending_loc = peptide[1]
            peptide_location = "[" + str(starting_loc) + " to " + str(ending_loc) + "]"
            SORTED_peptide_dict.setdefault(peptide_location, 0)
            SORTED_peptide_dict[peptide_location] = unique_peptide_list[counter]
            counter += 1

        return SORTED_peptide_dict



def BarGrapher(protein_name, oxidized_list, seq= None):
    if isinstance(oxidized_list, dict):
        plt.bar(range(len(oxidized_list)), list(oxidized_list.values()))
        plt.xticks(range(len(oxidized_list)), list(oxidized_list.keys()))
    else:
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

def write_to_excel(protein_name, oxidized_amino_acid_dict):

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
    print(file_name)
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
    '''Analysis'''
    #protein_name = "Polyubiquitin-B"
    protein_name = "Cytochrome c"
    Analyzer = PercentOxidationAnalyzer("CytC_2min.xlsx", 0, FASTA="ubiquitin.FASTA")
    Analyzer.PeptideValidator(protein_name)

    '''Save To Excel'''
    #write_to_excel(protein_name, oxidized_aminoacids)

if __name__ == "__main__":
    main()
