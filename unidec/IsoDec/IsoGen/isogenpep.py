import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogen_tools import peptide_to_dist, peptide_to_vector
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *


class IsoGenPepEngine(IsoGenEngineBase):
    def __init__(self, isolen=64):
        super().__init__()
        modelid=0
        self.isolen = isolen
        self.seqlengthranges = np.array([[1, 50], [51,200]])
        self.lengths = np.array([32, 64])
        self.inputname = "seqs"
        self.models = self.create_models()
        self.model = self.models[0]


    def create_models(self):
        #First create the peptide model
        pep = IsoGenModelBase(working_dir=None, isolen=64, vectorlen=31, savename="isogenpep_model.pth", modelid=0)
        pep.model = IsoGenNeuralNetwork(64,31)
        pep.savepath = os.path.join(pep.working_dir, pep.savename)
        pep.load_model()

        prot = IsoGenModelBase(working_dir=None, isolen=64, vectorlen=31, savename="isogenprot_model.pth", modelid=1)
        prot.model = IsoGenNeuralNetwork(64,31)
        prot.savepath = os.path.join(prot.working_dir, prot.savename)
        prot.load_model()

        return [pep, prot]


    def inputs_to_vectors(self, inputs):
        return np.array([peptide_to_vector(s) for s in inputs])


    def get_model_index(self, seqlength):
        if seqlength < 50:
            index = 0
        elif seqlength <= 200:
            index = 1
        else:
            print("Sequence length out of range")
            raise ValueError("Sequence length out of range, must be under 200 AAs.")
        return index


    def predict(self, seq):
        self.check(seq)
        vec = peptide_to_vector(seq)
        modelindex = self.get_model_index(len(seq))
        model = self.models[modelindex]
        return model.predict(vec)


    def check(self, seq):

        if len(seq) > 200:
            print("Sequence length longer than training data. Behavior may be unpredictable.")




if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    os.chdir("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Training")

    isolen = 32
    # trainfile = "peptidedists_633886.npz"
    trainfile = "Z:\\Group Share\\JGP\\PeptideTraining\\peptidedists_2492495.npz"
    # trainfile_synthetic = "peptidedists_synthetic_168420.npz"
    trainfile_synthetic = "Z:\\Group Share\\JGP\\PeptideTraining\\peptidedists_synthetic_3368720.npz"

    intact1 = "training_random_human_proteins_10000_min_50_max_200.npz"
    intact2 = "training_random_human_proteins_10000_min_50_max_200_mods2.npz"
    intact3 = "training_random_mouse_proteins_10000_min_50_max_200.npz"
    intact4 = "training_random_mouse_proteins_10000_min_50_max_200_mods2.npz"
    intact5 = "training_random_ecoli_proteins_10000_min_50_max_200.npz"
    intact6 = "training_random_ecoli_proteins_10000_min_50_max_200_mods2.npz"
    intact7 = "training_random_yeast_proteins_10000_min_50_max_200.npz"
    intact8 = "training_random_yeast_proteins_10000_min_50_max_200_mods2.npz"

    pep1 = "training_random_ecoli_proteins_10000_min_5_max_50_mods0.npz"
    pep2 = "training_random_ecoli_proteins_10000_min_5_max_50_mods1.npz"
    pep3 = "training_random_human_proteins_10000_min_5_max_50_mods0.npz"
    pep4 = "training_random_human_proteins_10000_min_5_max_50_mods1.npz"
    pep5 = "training_random_mouse_proteins_10000_min_5_max_50_mods0.npz"
    pep6 = "training_random_mouse_proteins_10000_min_5_max_50_mods1.npz"
    pep7 = "training_random_yeast_proteins_10000_min_5_max_50_mods0.npz"
    pep8 = "training_random_yeast_proteins_10000_min_5_max_50_mods1.npz"

    load = "Z:\\Group Share\\JGP\\PeptideTraining\\peptidedists_synthetic_3368720_short.npz"
    if False:
        print("Training Peptide Model...")
        eng_pep.train_multiple([trainfile_synthetic, pep1, pep2, pep3, pep4], epochs=10, forcenew=True)
        # print("Training Protein Model...")
        # eng_prot.train_multiple([intact1, intact2, intact3, intact4], epochs=10, forcenew=True)

    if False:
        testformulas = ["PEPTIDE", "CCCCCCCCCCCCC", "APTIGGGQGAAAAAAAAAAAASVGGTIPGPGPGGGQGPGEGGEGQTAR", "LLL", "KKK", "CCCM"]
        mpl.use("WxAgg")

        for i, f in enumerate(testformulas):
            maxval = 10
            dist = eng.predict(f)
            dist = dist / np.max(dist)
            truedist = peptide_to_dist(f)
            truedist = truedist/np.max(truedist)

            plt.subplot(2, 3, i + 1)
            plt.plot(dist * 100, label="AI", color="b")
            plt.plot(truedist * 100, color="k", label="True")
            plt.xlim(0, maxval)
            title_string = str(f)
            if len(f) > 8:
                title_string = title_string[:8] + "...{" + str(len(f)) + "}"
            plt.title(title_string)
            plt.xlabel("Isotope Number")
            plt.ylabel("%")
            if i == 3:
                plt.legend()
        plt.tight_layout()
        plt.show()

    if True:
        eng = IsoGenPepEngine()
        data = np.load("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Training\\human_protein_seqs.npz")

        dists = np.array(data["dists"])
        vecs = np.array(data["vecs"])
        seqs = np.array(data["seqs"])


        chisquareds = []

        for i in range(len(dists)):
            pred_dist = eng.predict(seqs[i])
            truedist = peptide_to_dist(seqs[i])
            #Truncate dists to the length of the predicted distribution
            chi = calculate_pearson_chisquared(pred_dist, dists[i])

            if i % 100 == 0:
                plt.plot(pred_dist, label="AI", color="b")
                plt.plot(truedist, color="r", label="True")
                plt.legend()
                plt.title(str(len(seqs[i])))
                plt.show()


            chisquareds.append(chi)

        chisquareds = np.array(chisquareds)
        print("Mean Chi Squared Error:", np.mean(chisquareds))
