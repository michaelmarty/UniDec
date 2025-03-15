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
        self.model = IsoGenModelBase(working_dir=None, isolen=isolen, vectorlen=20, savename="isogenpep_model_")
        self.isolen = isolen
        self.inputname = "seqs"

    def inputs_to_vectors(self, inputs):
        return np.array([peptide_to_vector(s) for s in inputs])

    def predict(self, seq):
        self.check(seq)
        vec = peptide_to_vector(seq)
        return self.model.predict(vec)

    def check(self, seq):
        if len(seq) > 300:
            print("Sequence length longer than training data. Behavior may be unpredictable.")




if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    os.chdir("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Training")

    isolen = 32
    trainfile = "peptidedists_633886.npz"
    trainfile = "peptidedists_2492495.npz"
    trainfile_synthetic = "peptidedists_synthetic_168420.npz"
    trainfile_synthetic = "Z:\\Group Share\\JGP\\PeptideTraining\\peptidedists_synthetic_3368720.npz"

    trainfile1 = "human_protein_seqs.npz"
    trainfile2 = "rat_protein_seqs.npz"
    trainfile3 = "yeast_protein_seqs.npz"
    trainfile4 = "b_subtilis_protein_seqs.npz"
    trainfile5 = "d_melanogaster_protein_seqs.npz"
    trainfile6 = "random_prots_yeast.npz"
    trainfile7 = "random_prots_human.npz"

    eng = IsoGenPepEngine(isolen=isolen)
    if True:
        eng.train_multiple([trainfile1,trainfile2,trainfile3,trainfile4,trainfile5, trainfile6, trainfile7],
                           epochs=10, forcenew=False)

    if True:
        peptide = peptide_to_vector("APTIGGGQGAAAAAAAAAAAASVGGTIPGPGPGGGQGPGEGGEGQTAR")

        print("Peptide:", peptide)
        exit()


        testformulas = ["PEPTIDE", "CCCCCCCCCCCCC", "APTIGGGQGAAAAAAAAAAAASVGGTIPGPGPGGGQGPGEGGEGQTAR", "LLL", "KKK", "CCCM"]
        mpl.use("WxAgg")
        for i, f in enumerate(testformulas):
            maxval = 10
            dist = eng.predict(f)
            truedist = peptide_to_dist(f)

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
