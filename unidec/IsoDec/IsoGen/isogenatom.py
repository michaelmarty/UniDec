import matplotlib.pyplot as plt
import matplotlib as mpl
from isogen_base import *
from isogen_tools import formula_to_vector, get_dist_from_formula


class IsoGenAtomEngine(IsoGenEngineBase):
    def __init__(self, isolen=128):
        super().__init__()
        self.model = IsoGenModelBase(working_dir=None, isolen=isolen, savename="isogenatom_model_", vectorlen=109)
        self.isolen = isolen
        self.inputname = "formulas"

    def inputs_to_vectors(self, inputs):
        return np.array([formula_to_vector(s) for s in inputs])

    def predict(self, formula):
        vec = formula_to_vector(formula)
        return self.model.predict(vec)


if __name__ == "__main__":
    os.chdir("Z:\\Group Share\\JGP\\PubChem")
    isolen = 32
    trainfile = "isodists_7056.npz"
    trainfile = "isodists_63202.npz"
    trainfile = "isodists_867971.npz"
    trainfile = "isodists_4194533.npz"

    trainfile_synethetic = "isodists_synthetic_1308425.npz"

    eng = IsoGenAtomEngine(isolen=isolen)
    if False:
        eng.train_multiple([trainfile, trainfile_synethetic], epochs=10, forcenew=False, inputname="formulas")

    testformulas = ["C2H5O", "H8N2O4W", "S8", "CHCl3", "C60", "CHBr3"]
    mpl.use("WxAgg")

    for i, f in enumerate(testformulas):
        maxval = 10
        dist = eng.predict(f)
        truedist = get_dist_from_formula(f, isolen=isolen)

        plt.subplot(2, 3, i + 1)
        plt.plot(dist * 100, label="AI", color="b")
        plt.plot(truedist * 100, color="k", label="True")
        plt.xlim(0, maxval)
        plt.title(str(f))
        plt.xlabel("Isotope Number")
        plt.ylabel("%")
        if i == 3:
            plt.legend()
    plt.tight_layout()
    plt.show()
