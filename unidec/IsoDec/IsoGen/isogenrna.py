import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenc import fft_rna_mass_to_dist, rna_to_vector, rna_vec_to_mass
from unidec.modules.isotopetools import isomike


class IsoGenRNAEngine(IsoGenEngineBase):
    def __init__(self, isolen=32):
        super().__init__()
        self.model = IsoGenModelBase(isolen=isolen, savename="isogenrna_model_", vectorlen=4)
        self.isolen = isolen
        self.raw_sequences = []

    def inputs_to_vectors(self, inputs):
        return np.array([rnaseq_to_vector(m) for m in inputs])

    def predict(self, seq):
        self.check(seq)
        vec = rnaseq_to_vector(seq)
        return self.model.predict(vec)

    def check(self, seq):
        if len(seq) > 100:
            print("Warning: Sequence is too long. Behavior may be unpredictable.")

    def scrape_rna_data(self, path):
        #masses = []
        dists = []
        seqs = []
        file = open(path, "r")

        for line in file:
            seq = line.strip()
            if seq[0] == '>':
                continue
            #self.raw_sequences.append(seq)
            vec = rna_to_vector(seq)
            mass = rna_vec_to_mass(vec)
            #masses.append(mass)
            seqs.append(seq)
            dists.append(fft_rna_mass_to_dist(mass))
        self.raw_sequences = np.array(self.raw_sequences)
        #np.savetxt("rna_raw_sequences.txt", self.raw_sequences, fmt = "%s", delimiter = "\n")
        return np.array(seqs), np.array(dists)




    def save_scraped_data(self, file_path, npz_outfile="rnadists.npz"):
        """
        Scrape RNA sequences, compute masses and isotopic distributions, then save as a .npz file.

        :param fasta_path: Path to the .fa file
        :param npz_outfile: Output path for the .npz dataset
        """
        self.outfile = npz_outfile
        X, Y = self.scrape_rna_data(file_path)
        print("Type of first X element:", type(X[0]))
        print("Example X[0]:", X[0])
        np.savez(npz_outfile, X=X, Y=Y)
        return X,Y

    def train_and_save_model(self, isolen, vectorlen, train_fname, epochs, forcenew=False):

        model = IsoGenModelBase(isolen=isolen, vectorlen=vectorlen, savename="isogenrna_model_")
        model.setup_model(forcenew=forcenew)
        model.setup_training()

        data = np.load(train_fname)
        keys = data.files
        if "masses" in keys and "dists" in keys:
            X, Y = data["masses"], data["dists"]
        elif "X" in keys and "Y" in keys:
            X, Y = data["X"], data["Y"]
        else:
            raise KeyError(f"Unexpected keys in {train_fname}: {keys}")

        train_loader, test_loader = self.setup_data(Y, X)

        model.run_training(train_loader, test_loader, epochs=epochs, forcenew=forcenew)
        model.save_model()



if __name__ == "__main__":
    path = "C:\\Users\\MartyLabsOfficePC\\Downloads\\ncbi_dataset\\ncbi_dataset\\data\\GCF_000001405.40\\rna.fna"
    eng = IsoGenRNAEngine(isolen=32)
    eng.save_scraped_data(path, npz_outfile="rna seqs.npz")


    #data = np.load("C:\\Python\\UniDec3\\unidec\\IsoDec\\IsoGen\\rnadists.npz")
    #Show the keys (like "X", "Y")
    # print("Keys:", data.files)
    #
    # # Access values
    # X = data["X"]
    # Y = data["Y"]
    # print("First 5 sequences (X):", X[:5])
    # print("First 5 distributions (Y):", Y[:5])
    exit()
    # Show first 5 entries of each

    isolens = [8,32]
    # for isolen in isolens:
    #     eng.train_and_save_model(isolen, train_fname="C:\\Python\\UniDec3\\unidec\\IsoDec\\IsoGen\\ncbi rna dataset.npz", vectorlen=4, epochs=10, forcenew=False)
    #
    # exit()

    os.chdir("Z:\\Group Share\\JGP\\PeptideTraining")

    for isolen in isolens:
        eng = IsoGenRNAEngine(isolen=isolen)

        trainfile_synthetic = "rnadists_synthetic_101424.npz"
        trainfile_synthetic2 =  "rna seqs.npz"


        if True:
            eng.train_multiple([trainfile_synthetic], epochs=10, forcenew=True)

    exit()
    mpl.use("WxAgg")
    # exit()
    testseqs = ["GUAC", "GUACGUACGUACGUAC", "GUACGUACGUACGUACGUACGUACGUACGUACGUAC", "GUACGUACGUACGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC", "AAAAAAAAAAAAA", "UUUUUUUUUUUUU"]
    #testmasses = [110000, 230000, 450000, 670000, 890000, 1000000]
    for i, m in enumerate(testseqs):
        dist = eng.predict(m)
        truedist = rnaseq_to_dist(m, isolen=isolen)
        plt.subplot(2, 3, i + 1)
        plt.plot(dist * 100, label="AI", color="b")
        plt.plot(truedist * 100, color="k", label="True")
        plt.xlim(0, isolen)
        title_string = str(m)
        if len(m) > 5:
            title_string = title_string[:5] + "...{" + str(len(m)) + "}"
        plt.title(title_string)
        plt.xlabel("Isotope Number")
        plt.ylabel("%")
        if i == 5:
            plt.legend()
    plt.tight_layout()
    plt.show()

