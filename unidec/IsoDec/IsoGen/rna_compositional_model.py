import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from unidec.IsoDec.IsoGen.isogen_tools import *


"""
Uses the compositional model to predict the isotopic distribution of a nucleic acid
The model is based on the work of D. Valkenborg from the paper "A Compositional Model to Predict the Aggregated Isotope Distribution
for Average DNA and RNA Oligonucleotides" metabolites, 2021

Limited to about 25 kDa!!
"""
class NucleicAcidCompositionalModelEngine:
    def __init__(self, analyte="RNA", isolen=20):
        self.isolen = isolen
        self.analyte = analyte
        self.k = 10 #number of polynomial components

        if analyte == "RNA":
            self.mu = 23151.52
            self.sigma = 4883.142
        elif analyte == "DNA":
            self.mu = 22746.2
            self.sigma = 4788.878
        else:
            raise ValueError("Analyte must be RNA or DNA")

        self.coefs = self.load_and_order_coefs(analyte=analyte)

        self.mass_shifts = np.array([0,1.002699,2.005361,3.007994,4.010602,
               5.013188,6.015755,7.018305,8.020840,
               9.023362,10.025870,11.028368,12.030855,
               13.033332,14.035800,15.038260,16.040712,
               17.043157,18.045595,19.048027])


    def load_and_order_coefs(self, analyte="RNA"):
        if analyte == "RNA":
            coefs = np.loadtxt("C:\\Python\\UniDec3\\unidec\\IsoDec\\IsoGen\\comp_model_parameters\\rna_polycoefs.txt")
        elif analyte == "DNA":
            coefs = np.loadtxt("C:\\Python\\UniDec3\\unidec\\IsoDec\\IsoGen\\comp_model_parameters\\dna_polycoefs.txt")

        ordered_coefs = [[0 for _ in range(20)] for _ in range(11)]

        current_col = 0
        current_row = 0
        total_vals = 0

        for i in range(len(coefs)):
            for j in range(len(coefs[i])):
                ordered_coefs[current_col][current_row] = coefs[i][j]
                total_vals += 1
                current_col += 1
                if current_col == 11:
                    total_vals = 0
                    current_col = 0
                    current_row += 1

        return ordered_coefs

    # Get the isotopic distribution from a mass
    def get_dist_from_mass(self, mass, charge=1, adductmass=-1.007276467):
        z = (mass - 23151.52) / 4883.142
        mass_mat = np.array([z for i in range(11)])
        for i in range(self.k+1):
            mass_mat[i] = mass_mat[i] ** i

        mass_mat = mass_mat.T

        q_ALR_hat = np.array(np.dot(mass_mat, self.coefs))
        q_ALR_hat = q_ALR_hat.reshape(-1, 1)  # Reshape to make it a 2D array (20, 1)
        print("shape of q_ALR_hat", q_ALR_hat.shape)
        # perform anti-ALR via softmax with one added to vector
        # Add a 1 to the beginning of q_ALR_hat and remove the last element
        tmp = np.insert(np.exp(q_ALR_hat), 0, 1, axis=0)
        tmp = np.delete(tmp, -1, axis=0)  # remove last element
        rs = np.sum(tmp, axis=0)  # row sums
        q_hat = tmp / rs

        #now create an m/z distribution
        masses = mass + self.mass_shifts
        mzs = np.abs((masses + adductmass) / charge)

        #combine the m/z and intensity arrays
        q_hat = np.array(q_hat).flatten()
        mzs = np.array(mzs).flatten()
        q_hat = q_hat / np.sum(q_hat)
        #Zip the two arrays together
        mzs = np.array(mzs)
        q_hat = np.array(q_hat)
        dist = np.column_stack((mzs, q_hat))
        return dist

    def get_intensities_from_mass(self, mass):
        z = (mass - 23151.52) / 4883.142
        mass_mat = np.array([z for i in range(11)])
        for i in range(self.k + 1):
            mass_mat[i] = mass_mat[i] ** i

        mass_mat = mass_mat.T

        q_ALR_hat = np.array(np.dot(mass_mat, self.coefs))
        q_ALR_hat = q_ALR_hat.reshape(-1, 1)  # Reshape to make it a 2D array (20, 1)
        # perform anti-ALR via softmax with one added to vector
        # Add a 1 to the beginning of q_ALR_hat and remove the last element
        tmp = np.insert(np.exp(q_ALR_hat), 0, 1, axis=0)
        tmp = np.delete(tmp, -1, axis=0)  # remove last element
        rs = np.sum(tmp, axis=0)  # row sums
        q_hat = tmp / rs

        return q_hat

if __name__ == "__main__":
    mpl.use("WxAgg")
    length = 500
    rna_seq_min = "U" * length
    rna_seq_max = "A" * length

    dist1 = rnaseq_to_dist(rna_seq_min, isolen=128)
    dist2 = rnaseq_to_dist(rna_seq_max, isolen=128)
    plt.plot(dist1, label="Min RNA")
    plt.plot(dist2, label="Max RNA")
    plt.show()





