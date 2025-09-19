import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenc import *


class IsoGenRNAEngine(IsoGenEngineBase):
    def __init__(self, isolen=64):
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
        if len(seq) > 300:
            print("Warning: Sequence is too long. Behavior may be unpredictable.")




if __name__ == "__main__":
    os.chdir("Z:\Group Share\JGP\PeptideTraining")
    trainfile = "rnadists_synthetic_1001424.npz"

    eng = IsoGenRNAEngine(isolen=64)
    #eng.train(trainfile, forcenew=True)

    mpl.use("WxAgg")

    testseqs = ["GUAC", "GUACGUACGUACGUAC", "GUACGUACGUACGUACGUACGUACGUACGUACGUAC", "GUACGUACGUACGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC", "AAAAAAAAAAAAA", "UUUUUUUUUUUUU"]
    #testmasses = [110000, 230000, 450000, 670000, 890000, 1000000]
    for i, m in enumerate(testseqs):
        dist = eng.predict(m)
        truedist = rnaseq_to_dist(m, isolen=64)
        plt.subplot(2, 3, i + 1)
        plt.plot(dist * 100, label="AI", color="b")
        plt.plot(truedist * 100, color="k", label="True")
        plt.xlim(0, 64)
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

