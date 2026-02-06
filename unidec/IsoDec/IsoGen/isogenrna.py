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

    trainfile = "synthetic_RNAs_10621.npz"
    trainfile1 = "training_random_RNAs_10000_min_21_max_220.npz"
    trainfile2 = "training_random_RNAs_10000_min_180_max_520.npz"

    eng1 = IsoGenRNAEngine(isolen=64)
    eng2 = IsoGenRNAEngine(isolen=128)

    #eng1.train_multiple([trainfile, trainfile1], epochs=10, forcenew=True)
    eng2.train_multiple([trainfile, trainfile1], epochs=10, forcenew=True)
    eng2.train(trainfile2, epochs=10, forcenew=False)


