import torch
import torch.nn
import numpy as np
import os
import platform
from torch import nn
from torch.optim import lr_scheduler
import inspect
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from numpy import fft as fftpack

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

mpl.use("WxAgg")

massavgine = 111.1254
avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])
isotopes = np.array(
    [
        [[12.0, 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
        [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
        [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
        [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
        [
            [31.9720707, 95.02],
            [32.9714585, 0.75],
            [33.9678668, 4.21],
            [35.9760809, 0.02],
        ],
    ]
)


# @njit(fastmath=True)
def makemass(testmass):
    num = testmass / massavgine * avgine
    intnum = np.array([int(round(n)) for n in num])
    x = intnum * isotopes[:, 0, 0]
    minmassint = np.sum(x)
    formula = ""
    if intnum[0] != 0:
        formula = formula + "C" + str(intnum[0])
    if intnum[1] != 0:
        formula = formula + "H" + str(intnum[1])
    if intnum[2] != 0:
        formula = formula + "N" + str(intnum[2])
    if intnum[3] != 0:
        formula = formula + "O" + str(intnum[3])
    if intnum[4] != 0:
        formula = formula + "S" + str(intnum[4])
    return formula, minmassint, intnum


# @njit(fastmath=True)
def isojim(isolist, length=128):
    """Thanks to Jim Prell for Sketching this Code"""
    buffer = np.zeros(length)
    h = np.array([1, 0.00015, 0, 0])
    c = np.array([1, 0.011, 0, 0])
    n = np.array([1, 0.0037, 0, 0])
    o = np.array([1, 0.0004, 0.002, 0])
    s = np.array([1, 0.0079, 0.044, 0])
    h = np.append(h, buffer)
    c = np.append(c, buffer)
    n = np.append(n, buffer)
    o = np.append(o, buffer)
    s = np.append(s, buffer)

    dt = np.dtype(np.complex128)
    hft = fftpack.rfft(h).astype(dt)
    cft = fftpack.rfft(c).astype(dt)
    nft = fftpack.rfft(n).astype(dt)
    oft = fftpack.rfft(o).astype(dt)
    sft = fftpack.rfft(s).astype(dt)

    numc = isolist[0]
    numh = isolist[1]
    numn = isolist[2]
    numo = isolist[3]
    nums = isolist[4]

    allft = cft**numc * hft**numh * nft**numn * oft**numo * sft**nums

    # with nb.objmode(allift='float64[:]'):
    #    allift = fftpack.irfft(allft)
    allift = np.abs(fftpack.irfft(allft))
    # allift = np.abs(allift)
    allift = allift / np.amax(allift)
    return allift[:length]  # .astype(nb.float64)


class IsoGenDataset(torch.utils.data.Dataset):
    """
    Dataset class for IsoGen
    """

    def __init__(self, masses, dists):
        """
        Initialize the dataset
        :param emat: List of encoded matrices
        :param z: List of charge assignments
        """
        self.dists = [torch.as_tensor(d, dtype=torch.float32) for d in dists]
        self.masses = torch.as_tensor(masses, dtype=torch.float32)

    def __len__(self):
        return len(self.masses)

    def __getitem__(self, idx):
        return [self.masses[idx], self.dists[idx]]


def true_dist(mass, isolength=128):
    _, minmassint, isolist = makemass(mass)
    intensities = isojim(isolist, length=isolength)
    return intensities


def gen_training_data(n, isolen=128, massrange=[100, 1000000], log=False):
    if log:
        rlog = np.random.uniform(np.log10(massrange[0]), np.log10(massrange[1]), n)
        randmasses = np.power(10, rlog)
    else:
        randmasses = np.random.uniform(massrange[0], massrange[1], n)
    dists = [true_dist(m, isolen) for m in randmasses]
    return randmasses, np.array(dists)


def create_data_loaders(training_data, test_data):
    training_data = IsoGenDataset(training_data[0], training_data[1])
    test_data = IsoGenDataset(test_data[0], test_data[1])

    train_dataloader = DataLoader(
        training_data, batch_size=32, shuffle=True, pin_memory=True
    )
    test_dataloader = DataLoader(
        test_data, batch_size=2048, shuffle=False, pin_memory=False
    )

    return train_dataloader, test_dataloader


class IsoGenModel:
    """
    General model class for isotope distrubution generation.

    Includes functions to train, evaluate, and predict isotope distributions.
    """

    def __init__(self, working_dir=None, isolen=64):
        """
        Initialize the model. Set the working directory
        :param working_dir:
        """
        if working_dir is None:
            if platform.system() == "Linux":
                # working_dir = "/xdisk/mtmarty/mtmarty/training/"
                filename = inspect.getframeinfo(inspect.currentframe()).filename
                working_dir = os.path.dirname(os.path.abspath(filename))
            else:
                # working_dir = "C:\\Data\\IsoNN\\"
                filename = inspect.getframeinfo(inspect.currentframe()).filename
                working_dir = os.path.dirname(os.path.abspath(filename))

        self.working_dir = working_dir

        self.device = None
        self.model = None
        self.loss_fn = None
        self.optimizer = None
        self.scheduler = None

        self.modelid = None
        self.savepath = None
        self.indexes = None
        self.isolen = isolen
        self.modelid = 0
        # self.get_model(self.modelid)

    def get_model(self, modelid):
        """
        Get the model based on the model ID.
        :param modelid: Model ID integer. Options are 0, 1, and 2.
        :return: None
        """
        self.modelid = modelid
        if modelid == 0:
            self.model = IsoGenNeuralNetwork(isolen=self.isolen)
            savename = "isogen_model_" + str(self.isolen) + ".pth"
        else:
            print("Model ID not recognized", modelid)
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def setup_model(self, modelid=None, forcenew=False):
        """
        Setup model and load if savepath exists. Set device.
        :param modelid: Model ID passed to self.get_model()
        :param forcenew: Whether to force starting over from scratch on model parameters
        :return: None
        """
        if modelid is None:
            modelid = self.modelid
        self.device = (
            "cuda"
            if torch.cuda.is_available()
            else "mps" if torch.backends.mps.is_available() else "cpu"
        )
        print(f"Using {self.device} device")

        self.get_model(modelid)
        if os.path.isfile(self.savepath) and not forcenew:
            self.load_model()
            print("Loaded Weights:", self.savepath)
        else:
            print("Starting Model From Scratch. Let's ride!")
        self.model = self.model.to(self.device)
        print("Loaded Model:", self.model)

    def setup_training(self, lossfn="mse", forcenew=False):
        """ "
        Setup loss function, optimizer, and scheduler.
        :param lossfn: Loss function to use. Options are "crossentropy", "weightedcrossentropy", and "focal".
        :return: None
        """
        if self.model is None:
            self.setup_model(forcenew=forcenew)
        if lossfn == "mse":
            self.loss_fn = nn.MSELoss()
        else:
            raise ValueError("Loss function not recognized.", lossfn)
        # self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        # self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)
        self.scheduler = None

    def load_model(self):
        """
        Load model from savepath.
        :return: None
        """
        if os.path.isfile(self.savepath):
            try:
                self.model.load_state_dict(torch.load(self.savepath, weights_only=True))
                print("Model loaded:", self.savepath)
                # print_model(self.model)
            except Exception as e:
                print("Model failed to load:", self.savepath)
                print(e)
        else:
            print("Model not found:", self.savepath)

    def save_model(self):
        """
        Save model to savepath.
        :return: None
        """
        if self.model is None:
            self.setup_model()
        torch.save(self.model.state_dict(), self.savepath)
        # save_model_to_binary(self.model, self.savepath.replace(".pth", ".bin"))
        print("Model saved:", self.savepath, self.savepath.replace(".pth", ".bin"))

    def train_model(self, dataloader, lossfn="mse", forcenew=False):
        """
        Train the model on a DataLoader object.
        :param dataloader: Training DataLoader object
        :param lossfn: Loss function to use. Options are "crossentropy", "weightedcrossentropy", and "focal".
        :return: None
        """
        self.setup_training(lossfn, forcenew=forcenew)
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.train()

        for batch, s in enumerate(dataloader):
            x = s[0]
            y = s[1]
            x, y = x.to(self.device), y.to(self.device)
            x = x.unsqueeze(1)
            # Compute prediction error
            pred = self.model(x)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % int(num_batches / 5) == 0:
                loss, current = loss.item(), (batch + 1) * len(x)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
            if self.scheduler is not None:
                self.scheduler.step()

    def evaluate_model(self, dataloader):
        """
        Evaluate the model on a test set.
        :param dataloader: Test DataLoader object
        :return: None
        """
        if self.model is None:
            self.setup_training()

        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss = 0
        rmsdlist = []
        with torch.no_grad():
            for j, s in enumerate(dataloader):
                x = s[0]
                y = s[1]

                x, y = x.to(self.device), y.to(self.device)
                x = x.unsqueeze(1)
                pred = self.model(x)
                test_loss += self.loss_fn(pred, y).item()

                rmsd = (pred - y) ** 2
                rmsd = rmsd.cpu().detach().numpy()
                rmsd = np.sum(rmsd, axis=1)
                rmsdlist.append(rmsd)

        test_loss /= num_batches

        rmsdlist = np.concatenate(rmsdlist)
        avgrmsd = np.sum(rmsdlist) / len(rmsdlist)

        print(
            f"Test Error: \n Avg RMSD: {(avgrmsd*100):>0.1f}, Avg loss: {test_loss:>8f} \n"
        )

    def predict(self, mass):
        """
        Predict the isotope distribution for a given mass
        :param mass: Mass value
        :return: Predicted isotope distribution intensity vector
        """
        if self.model is None:
            self.setup_model()

        test_data = torch.tensor([mass], dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predvec = predvec.squeeze()

        # print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return predvec.cpu().detach().numpy()

    def batch_predict(self, masses):
        """
        Predict the isotope distribution for a list of masses
        :param masses: List of mass values
        :return: Predicted isotope distribution intensity vectors
        """
        if self.model is None:
            self.setup_model()

        test_data = torch.tensor(masses, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(1)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predvec = predvec.squeeze()

        return predvec.cpu().detach().numpy()

    def run_training(
        self,
        train_dataloader,
        test_dataloader,
        epochs=10,
        lossfn="mse",
        forcenew=False,
        save=True,
    ):
        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.train_model(train_dataloader, lossfn=lossfn, forcenew=forcenew)
            self.evaluate_model(test_dataloader)

        if save:
            self.save_model()


class IsoGenNeuralNetwork(nn.Module):
    """
    Very simple neural net for generating the isotope distirbution.
    """

    def __init__(self, isolen=128):
        super().__init__()
        h1 = isolen
        h2 = isolen
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(5, h1),
            nn.Tanh(),
            nn.Linear(h1, h2),
            nn.Tanh(),
            nn.Linear(h2, isolen),
            nn.Sigmoid(),
        )

    def forward(self, x):
        newx = self.encode(x)
        logits = self.linear_relu_stack(newx)
        return logits

    def encode(self, x):
        x1 = x / 1000000.0
        x2 = x / 100000.0
        x3 = x / 10000.0
        x4 = x / 1000.0
        x5 = x / 100.0
        newx = torch.concatenate((x1, x2, x3, x4, x5), dim=1)
        newx = torch.clamp(newx, 0, 1)
        return newx


class IsoGenEngine:
    def __init__(self):
        self.massranges = np.array(
            [[10, 1200], [800, 12000], [8000, 120000], [80000, 1000000]]
        )
        self.lengths = np.array([8, 32, 128, 1024])

        self.models = []
        for l in self.lengths:
            model = IsoGenModel(isolen=l)
            self.models.append(model)

    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
        except:
            print("Length not found in list")
            return
        massrange = self.massranges[modelindex]
        model = self.models[modelindex]
        input, target = gen_training_data(n, isolen=length, massrange=massrange)
        testinput, testtarget = gen_training_data(
            int(n * 0.1), isolen=length, massrange=massrange
        )
        print("Created Data:", length, n)
        trd, ted = create_data_loaders([input, target], [testinput, testtarget])
        model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

    def train_all(self, n=100000, epochs=10, forcenew=False):
        for i, l in enumerate(self.lengths):
            self.train(n, epochs, l)

    def get_isolen(self, mass):
        if mass < 1000:
            index = 0
        elif mass <= 10000:
            index = 1
        elif mass <= 100000:
            index = 2
        elif mass <= 1000000:
            index = 3
        else:
            print("Mass out of range", mass)
            raise ValueError("Mass out of range. Must be less than 1 MDa.")

        isolen = self.lengths[index]
        return isolen, index

    def predict(self, mass):
        isolen, index = self.get_isolen(mass)
        model = self.models[index]
        return model.predict(mass)


if __name__ == "__main__":

    eng = IsoGenEngine()

    if True:
        n = 100000
        #eng.train_all(n, 10, forcenew=False)
        #eng.train(n, 10, 128, forcenew=False)
        eng.train(n, 10, 1024, forcenew=False)


    # exit()
    testmasses = [100, 1000, 10000, 50000, 100000, 1000000]
    testmasses = [110000, 230000, 450000, 670000, 890000, 1000000]
    for i, m in enumerate(testmasses):
        isolen = eng.get_isolen(m)[0]
        maxval = isolen
        dist = eng.predict(m)
        truedist = true_dist(m, isolength=isolen)

        plt.subplot(2, 3, i + 1)
        plt.plot(dist * 100, label="AI", color="b")
        plt.plot(truedist * 100, color="k", label="True")
        plt.xlim(0, maxval)
        plt.title(str(m) + " Da")
        plt.xlabel("Isotope Number")
        plt.ylabel("%")
        if i == 3:
            plt.legend()
    plt.tight_layout()
    plt.show()
