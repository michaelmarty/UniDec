import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogen_tools import peptide_to_dist, peptide_to_vector
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *


class IsoGenPepEngine(IsoGenEngineBase):
    def __init__(self, isolen=32):
        super().__init__()
        modelid=0
        self.model = IsoGenModelBase(working_dir=None, isolen=isolen, vectorlen=31, savename="isogenpep_model_", modelid=modelid)
        self.isolen = isolen
        self.seqlengthranges = np.array([[1, 50], [50,200]])
        self.lengths = np.array([32, 64])
        self.inputname = "seqs"

        self.models = []
        for l in self.lengths:
            model = IsoGenModelBase(working_dir=None, isolen=l, vectorlen=20, savename="isogenpep_model_", modelid=modelid)
            self.models.append(model)


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
    intact2 = "training_random_yeast_proteins_10000_min_50_max_200.npz"
    intact3 = "training_random_mouse_proteins_10000_min_50_max_200.npz"
    intact4 = "training_random_ecoli_proteins_10000_min_50_max_200.npz"

    eng_pep = IsoGenPepEngine(isolen=32)
    eng_prot = IsoGenPepEngine(isolen=64)
    if True:
        # print("Training Peptide Model")
        # eng_pep.train_multiple([trainfile, trainfile_synthetic], epochs=10, forcenew=True)
        print("Training Protein Model")
        eng_prot.train_multiple([intact1, intact2, intact3, intact4], epochs=10, forcenew=True)

    eng = IsoGenPepEngine()
    if True:
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

    if False:
        data = np.load("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Evaluation\\ecoli_protein_seqs.npz")


        dists = np.array(data["dists"])
        vecs = np.array(data["vecs"])
        seqs = np.array(data["seqs"])


        data1 = np.load("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Evaluation\\mouse_protein_seqs.npz")

        dists = np.concatenate((dists, np.array(data1["dists"])))
        vecs = np.concatenate((vecs, np.array(data1["vecs"])))
        seqs = np.concatenate((seqs, np.array(data1["seqs"])))

        chisquareds = []

        for i in range(len(dists)):
            pred_dist = eng.predict(seqs[i])
            #Truncate dists to the length of the predicted distribution
            chi = calculate_pearson_chisquared(pred_dist, dists[i])

            if i % 100 == 0:
                plt.plot(pred_dist, label="AI", color="b")
                plt.plot(dists[i], color="r", label="True")
                plt.legend()
                plt.title(str(len(seqs[i])))
                plt.show()


            chisquareds.append(chi)

        chisquareds = np.array(chisquareds)
        print("Mean Chi Squared Error:", np.mean(chisquareds))
