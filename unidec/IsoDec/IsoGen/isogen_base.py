import torch
import numpy as np
import os
import platform
from torch import nn
from torch.optim import lr_scheduler
import inspect
from torch.utils.data import DataLoader
import time



os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

# Function to save the model to a binary format that can be read into C
def save_model_to_binary(model, outfile):
    params = model.parameters()
    output = []
    for m in params:
        flat = m.flatten()
        output.append(flat)

    output = torch.cat(output).cpu().detach().numpy()

    output.tofile(outfile)

class IsoGenDatasetBase(torch.utils.data.Dataset):
    """
    Dataset class for IsoGen
    """
    def __init__(self, inputs, dists, vectors):
        """
        Initialize the dataset
        :param emat: List of encoded matrices
        :param z: List of charge assignments
        """
        self.dists = [torch.as_tensor(d, dtype=torch.float32) for d in dists]
        self.inputs = inputs
        self.vectors = [torch.as_tensor(v, dtype=torch.float32) for v in vectors]

    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, idx):
        return [self.vectors[idx], self.dists[idx]]


class IsoGenModelBase:
    """
    General model class for isotope distrubution generation.

    Includes functions to train, evaluate, and predict isotope distributions.
    """

    def __init__(self, working_dir=None, isolen=128, vectorlen=20, savename="isogenpep_model_", modelid=0):
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
        self.savename = savename
        self.indexes = None
        self.isolen = isolen
        self.vectorlen = vectorlen
        self.modelid = modelid

    def get_model(self, modelid):
        """
        Get the model based on the model ID.
        :param modelid: Model ID integer. Options are 0, 1, and 2.
        :return: None
        """
        self.modelid = modelid
        if modelid == 0:
            self.model = IsoGenNeuralNetwork(isolen=self.isolen, vectorlen=self.vectorlen)
            savename = self.savename + str(self.isolen) + ".pth"
        elif modelid == 1:
            self.model = HighMassNeuralNetwork(isolen=self.isolen, vectorlen=self.vectorlen)
            savename = self.savename + str(self.isolen) + ".pth"
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
            self.loss_fn = nn.MSELoss(reduction="sum")
        else:

            raise ValueError("Loss function not recognized.", lossfn)

        # self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        #self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)
        self.scheduler = None

    def load_model(self):
        """
        Load model from savepath.
        :return: None
        """
        if os.path.isfile(self.savepath):
            try:
                if torch.cuda.is_available():
                    self.model.load_state_dict(torch.load(self.savepath, weights_only=True))
                    print("Model loaded:", self.savepath)
                    # print_model(self.model)
                else:
                    self.model.load_state_dict(torch.load(self.savepath, weights_only=True, map_location=torch.device('cpu')))
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
        save_model_to_binary(self.model, self.savepath.replace(".pth", ".bin"))
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

        log_interval = max(1, int(num_batches / 5))
        self.model.train()

        for batch, s in enumerate(dataloader):
            x = s[0]
            y = s[1]
            y = y.squeeze(1)
            x, y = x.to(self.device), y.to(self.device)
            x = x.unsqueeze(1)
            # Compute prediction error
            pred = self.model(x).squeeze(1)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % log_interval == 0:
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
                pred = self.model(x).squeeze(1)
                test_loss += self.loss_fn(pred, y).item()

                rmsd = (pred - y) ** 2
                rmsd = rmsd.cpu().detach().numpy()
                rmsd = np.sum(rmsd, axis=1)
                rmsdlist.append(rmsd)

        test_loss /= num_batches

        rmsdlist = np.concatenate(rmsdlist)
        avgrmsd = np.sum(rmsdlist) / len(rmsdlist)

        print(
            f"Test Error: \n Avg RMSD: {(avgrmsd*100):>0.2f}, Avg loss: {test_loss:>8f} \n"
        )

    def predict(self, vec):
        """
        Predict the isotope distribution for a given mass
        :param mass: Mass value
        :return: Predicted isotope distribution intensity vector
        """
        if self.model is None:
            self.setup_model()

        test_data = torch.tensor(vec, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predvec = predvec.squeeze()

        # print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return predvec.cpu().detach().numpy()

    def batch_predict(self, vectors):
        """
        Predict the isotope distribution for a list of masses
        :param masses: List of mass values
        :return: Predicted isotope distribution intensity vectors
        """
        if self.model is None:
            self.setup_model()

        test_data = torch.tensor(vectors, dtype=torch.float32)

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

    def __init__(self, isolen=128, vectorlen=20):
        super().__init__()
        h1 = isolen
        h2 = isolen
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(vectorlen, h1),
            nn.Softsign(),
            nn.Linear(h1, h2),
            nn.Softsign(),
            nn.Linear(h2, isolen),
            nn.Sigmoid(),
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits


class HighMassNeuralNetwork(nn.Module):
    """
    Very simple neural net for generating the isotope distirbution.
    """

    def __init__(self, isolen=1024, vectorlen=6):
        super().__init__()
        h1 = 32
        h2 = 256
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(vectorlen, h1),
            nn.PReLU(),
            nn.Linear(h1, h2),
            nn.PReLU(),
            nn.Linear(h2, isolen),
            nn.PReLU(),
            nn.Linear(isolen, isolen),
            nn.Tanh(),
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        #logits = logits / torch.amax(logits, dim=1, keepdim=True)
        return (logits+1)/2



class IsoGenEngineBase:
    def __init__(self, isolen=128):
        self.model = None
        self.isolen = isolen
        self.DatasetType = IsoGenDatasetBase
        self.inputname = "seqs"

    def inputs_to_vectors(self, inputs):
        return inputs

    def create_data_loaders(self, traindata, testdata, batchsize=32, testbatchsize=2048):
        vectors1 = self.inputs_to_vectors(traindata[0])
        vectors2 = self.inputs_to_vectors(testdata[0])
        print("Length of traindata[0]", len(traindata[0]))
        print("Length of traindata[1]", len(traindata[1]))
        training_data = self.DatasetType(traindata[0], traindata[1], vectors1)
        test_data = self.DatasetType(testdata[0], testdata[1], vectors2)

        train_dataloader = DataLoader(
            training_data, batch_size=batchsize, shuffle=True, pin_memory=True
        )
        test_dataloader = DataLoader(
            test_data, batch_size=testbatchsize, shuffle=False, pin_memory=False
        )
        print("Data Loaded:", len(training_data), len(test_data))
        return train_dataloader, test_dataloader

    def setup_data(self, dists, seqs, trainper = 0.9):
        if dists.shape[1] < self.isolen:
            print("Warning: Isolen does not match training data. Changing to", dists.shape[1])
            self.isolen = dists.shape[1]

        elif dists.shape[1] > self.isolen:
            print("Warning: Isolen does not match training data. Truncating to", self.isolen)
            dists = dists[:, :self.isolen]

        # Split randomly into training and testing with 90% training
        n = len(seqs)
        ntrain = int(n * trainper)
        indexes = np.arange(n)
        np.random.shuffle(indexes)
        trainindexes = indexes[:ntrain]
        testindexes = indexes[ntrain:]

        trainformulas = seqs[trainindexes]
        traindists = dists[trainindexes]

        testformulas = seqs[testindexes]
        testdists = dists[testindexes]

        trd, ted = self.create_data_loaders([trainformulas, traindists], [testformulas, testdists])
        return trd, ted

    def train(self, train_fname, epochs=10, forcenew=False, inputname=None):
        if inputname is None:
            inputname = self.inputname
        tdata = np.load(train_fname)
        dists = tdata["dists"]
        seqs = tdata[inputname]

        trd, ted = self.setup_data(dists, seqs)
        self.model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

    def train_multiple(self, train_fnames, epochs=10, forcenew=False, inputname=None):
        if inputname is None:
            inputname = self.inputname
        seqall = []
        distall = []
        for train_fname in train_fnames:
            tdata = np.load(train_fname)
            dists = tdata["dists"]
            seqs = tdata[inputname]
            seqall.extend(seqs)
            distall.extend(dists)

        trd, ted = self.setup_data(np.array(distall), np.array(seqall))
        self.model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

    def predict(self, seq):
        return None

    def check(self, seq):
        return None

if __name__ == "__main__":
    model = IsoGenNeuralNetwork(isolen=8, vectorlen=4)
    model.load_state_dict(torch.load("isogenrna_model_8.pth", map_location='cpu'))
    save_model_to_binary(model, "isogenrna_model_8.bin")


    # model = IsoGenNeuralNetwork(isolen=8, vectorlen=5)
    # model.load_state_dict(torch.load("isogenmass_model_8.pth", map_location='cpu'))
    # save_model_to_binary(model, "isogenmass_model_8.bin")
