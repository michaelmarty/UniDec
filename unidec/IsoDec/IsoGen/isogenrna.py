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
        self.isolen = isolen
        self.seqlengthranges = np.array([[1, 200], [201, 500]])
        self.lengths = np.array([64, 128])
        self.inputname = "seqs"
        self.models = []
        for l in self.lengths:
            if l > 128:
                modelid = 1
            else:
                modelid = 0
            model = IsoGenModelBase(isolen=l, savename="isogenrna_model_", vectorlen=4, modelid=modelid)
            self.models.append(model)

        for i in range(len(self.lengths)):
            if self.lengths[i] == self.isolen:
                self.model = self.models[i]

    def inputs_to_vectors(self, inputs):
        return np.array([rnaseq_to_vector(m) for m in inputs])

    def get_model_index(self, seqlength):
        for i, range in enumerate(self.seqlengthranges):
            if seqlength >= range[0] and seqlength <= range[1]:
                return i

        print("Sequence length out of range.")
        raise ValueError("Sequence length out of range, must be under 500.")

    def predict(self, seq):
        self.check(seq)
        vec = rnaseq_to_vector(seq)
        modelindex = self.get_model_index(len(seq))
        model = self.models[modelindex]
        return model.predict(vec)

    def check(self, seq):
        if len(seq) > 500:
            print("Warning: Sequence is too long. Behavior may be unpredictable.")




if __name__ == "__main__":
    os.chdir(r"C:\Users\Admin\Documents\martylab\RNA_SeqData\Training")

    trainfile0 = "synthetic_RNAs_min2_max35.npz"
    trainfile1 = "synthetic_RNAs_min36_max220.npz"
    trainfile2 = "synthetic_RNAs_min30_max120.npz"

    trainfile3 = "synthetic_RNAs_min180_max520.npz"
    trainfile4 = "synthetic_RNAs_min200_max400.npz"

    eng1 = IsoGenRNAEngine(isolen=64)
    eng2 = IsoGenRNAEngine(isolen=128)

    #eng1.train_multiple([trainfile0, trainfile0, trainfile1, trainfile2], epochs=20, forcenew=True)
    eng2.train_multiple([trainfile3, trainfile4], epochs=20, forcenew=True)