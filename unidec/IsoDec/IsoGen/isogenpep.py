import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from isogen_base import *
from isogen_tools import peptide_to_dist, peptide_to_vector


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
        if len(seq) > 40:
            print("Sequence length longer than training data. Behavior may be unpredictable.")


if __name__ == "__main__":
    os.chdir("Z:\\Group Share\\JGP\\PeptideTraining")

    isolen = 32
    trainfile = "peptidedists_633886.npz"
    trainfile = "peptidedists_2492495.npz"
    trainfile_synthetic = "peptidedists_synthetic_168420.npz"
    trainfile_synthetic = "peptidedists_synthetic_3368720.npz"

    eng = IsoGenPepEngine(isolen=isolen)
    if False:
        eng.train_multiple([trainfile, trainfile_synthetic], epochs=10, forcenew=False)

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
