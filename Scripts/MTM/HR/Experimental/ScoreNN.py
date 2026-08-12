import matplotlib.pyplot as plt

import time
import numpy as np
import os
import torch
from torch.utils.data import DataLoader
from torch import nn
from torch.optim import lr_scheduler
from score_encoding import vec_from_assignment, get_dummy_charge_state

import platform
import inspect

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


class ScoreNNModel:
    """
    General model class for isotope distrubution generation.

    Includes functions to train, evaluate, and predict isotope distributions.
    """

    def __init__(self, working_dir=None, vectorlen=32, savename="score_model_", modelid=0):
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
            self.model = ScoreNeuralNetwork(vectorlen=self.vectorlen)
            savename = self.savename + str(self.vectorlen) + ".pth"
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
        :param lossfn: Loss function to use. Options are "mse"
        :return: None
        """
        if self.model is None:
            self.setup_model(forcenew=forcenew)
        if lossfn == "mse":
            # self.loss_fn = nn.MSELoss(reduction="sum")
            self.loss_fn = nn.BCEWithLogitsLoss()
        else:
            raise ValueError("Loss function not recognized.", lossfn)
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        # self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.1)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)
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
        save_model_to_binary(self.model, self.savepath.replace(".pth", ".bin"))
        print("Model saved:", self.savepath, self.savepath.replace(".pth", ".bin"))

    def train_model(self, dataloader, lossfn="mse", forcenew=False):
        """
        Train the model on a DataLoader object.
        :param dataloader: Training DataLoader object
        :param lossfn: Loss function to use. Options are "mse".
        :return: None
        """
        self.setup_training(lossfn, forcenew=forcenew)
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.train()

        for batch, s in enumerate(dataloader):
            x = s[0]
            y = s[1]
            y = y.unsqueeze(1)
            x, y = x.to(self.device), y.to(self.device)
            x = x.unsqueeze(1)
            # Compute prediction error
            pred = self.model(x).squeeze(1)
            loss = self.loss_fn(pred, y)

            # Test if loss is nan
            if torch.isnan(loss):
                print(batch)
                print(x)
                print(y)
                print(pred)

                # Find the index in pred that is nan
                b1 = torch.isnan(pred)
                bad = x[b1]
                print(bad)

                raise ValueError("Loss is NaN")

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

        num_batches = len(dataloader)
        self.model.eval()
        test_loss = 0

        ncorrect = 0
        ntotal = 0
        with torch.no_grad():
            for j, s in enumerate(dataloader):
                x = s[0]
                y = s[1]

                x, y = x.to(self.device), y.to(self.device)
                x = x.unsqueeze(1)
                y = y.unsqueeze(1)
                pred = self.model(x).squeeze(1)
                test_loss += self.loss_fn(pred, y).item()

                # Calculate accuracy
                sigmoid = nn.Sigmoid()
                pred = sigmoid(pred)
                pred = pred > 0.5
                correct = pred == y
                ncorrect += correct.sum().item()
                ntotal += len(y)

        test_loss /= num_batches

        print(f"Test Error: \n Avg loss: {test_loss:>8f} Percent Correct: {ncorrect/ntotal * 100}% \n")

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
        sigmoid = nn.Sigmoid()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predvec = predvec.squeeze()
            predvec = sigmoid(predvec)

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
        sigmoid = nn.Sigmoid()

        x = test_data.unsqueeze(1)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predvec = predvec.squeeze()
            predvec = sigmoid(predvec)

        return predvec.cpu().detach().numpy()


class ScoreNeuralNetwork(nn.Module):
    """
    Very simple neural net for generating the isotope distirbution.
    """

    def __init__(self, vectorlen=32):
        super().__init__()
        h1 = 32
        h2 = 16
        h3 = 8
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(vectorlen, h1),
            nn.Softsign(),
            nn.Linear(h1, h2),
            nn.Softsign(),
            nn.Linear(h2, h3),
            nn.Softsign(),
            nn.Linear(h3, 1),
            # nn.Sigmoid(),
        )

    def forward(self, x):
        x = torch.abs(x)
        logits = self.linear_relu_stack(x)
        return logits


class SNNDataSet(torch.utils.data.Dataset):
    """
    Dataset class for ScoreNN
    """

    def __init__(self, vec, res):
        """
        Initialize the dataset
        """
        self.vec = [torch.as_tensor(v, dtype=torch.float32) for v in vec]
        self.res = torch.as_tensor(res, dtype=torch.float32)
        self.check_data()

    def __len__(self):
        return len(self.res)

    def __getitem__(self, idx):
        return [self.vec[idx], self.res[idx]]

    def check_data(self):
        print("Checking min and max of res:", torch.min(self.res), torch.max(self.res))
        badindexes = []
        for i, v in enumerate(self.vec):
            if torch.isnan(v).any():
                badindexes.append(i)
                # raise ValueError("NAN in vec")

        if len(badindexes) > 0:
            # Delete bad indexes
            print("Deleting bad indexes:", badindexes)
            self.vec = [v for i, v in enumerate(self.vec) if i not in badindexes]
            self.res = [r for i, r in enumerate(self.res) if i not in badindexes]


class ScoreNNEngine:
    """
    Main class for IsoDec Engine
    """

    def __init__(self):
        """
        Initialize the Engine
        """
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.test_centroids = []
        self.training_centroids = []
        self.scoremodel = ScoreNNModel()
        self.batch_size = 32
        self.test_batch_size = 1024
        self.isolength = 16

    def load_training_data(self, training_path, test_path=None):
        """
        Load training data from a file
        :param training_path: Path to the training data file or name of the file tag
        :param test_path: Optional path to the test data file or name of the file tag.
                            If not, will default to same as training_path
        :param noise_percent: The percent of noise to add to the training and test data
        :param double_percent: The percent of double peaks to add to the training and test data
        :return: None
        """
        ext = ".npz"
        if ext not in training_path:
            training_path = "training_data_" + training_path + ext

        if test_path is None:
            test_path = training_path.replace("training", "test")
        elif ext not in test_path:
            test_path = "test_data_" + test_path + ext

        td = np.load(training_path, allow_pickle=True)
        self.training_data = [td["vecs"], td["res"]]
        td = np.load(test_path, allow_pickle=True)
        self.test_data = [td["vecs"], td["res"]]

        print("Loaded:", len(self.training_data[0]), "Training Samples")

    def create_merged_dataloader(self, dirs, training_path, batchsize=None):
        """
        Create a merged dataloader from multiple directories. Looks for common file names and merges them together
        :param dirs: Directories to look in
        :param training_path: File name or tag, fed to load_training_data
        :param noise_percent: Percent of noise to add to the training and test data
        :param batchsize: Batch size for training
        :param double_percent: Percent of double peaks to add to the training and test data
        :return:
        """
        if batchsize is not None:
            self.batch_size = batchsize

        training_data = []
        test_data = []

        for d in dirs:
            os.chdir(d)
            print(d)
            self.load_training_data(training_path)
            training_data.append(self.training_data)
            test_data.append(self.test_data)

        self.training_data = [np.concatenate([t[0] for t in training_data], axis=0),
                              np.concatenate([t[1] for t in training_data], axis=0)]
        self.test_data = [np.concatenate([t[0] for t in test_data], axis=0),
                          np.concatenate([t[1] for t in test_data], axis=0)]

        self.training_data = SNNDataSet(self.training_data[0], self.training_data[1])
        self.test_data = SNNDataSet(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.test_batch_size, shuffle=False,
                                          pin_memory=True)

        print(f"Training Data Length: {len(self.training_data)}")
        print(f"Test Data Length: {len(self.test_data)}")

    def train_model(self, epochs=30, save=True, lossfn="mse", forcenew=False):
        """
        Train the model
        :param epochs: Number of epochs
        :param save: Whether to save it. Default is True
        :param lossfn: Loss function, default is mse
        :return: None
        """
        starttime = time.perf_counter()

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created.")
        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.scoremodel.train_model(self.train_dataloader, lossfn=lossfn, forcenew=forcenew)
            self.scoremodel.evaluate_model(self.test_dataloader)

        if save:
            self.scoremodel.save_model()
        print("Done! Time:", time.perf_counter() - starttime)

    def predictor(self, centroids, z):
        """
        Predict the score for a given charge state assignment
        :param centroids: Set of centroid data for a peak with m/z in first column and intensity in second
        :param z: Charge state assignment
        :return: Score
        """
        cv = vec_from_assignment(centroids, z, isolen=self.isolength)
        score = self.scoremodel.predict(cv)
        return score


if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = ScoreNNEngine()
    # topdirectory = "C:\\Data\\IsoNN\\training"
    #
    # dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]
    # eng.create_merged_dataloader(dirs, "phase82")

    eng.create_merged_dataloader(["Z:\\Group Share\\JGP\\MSV000090488\\"], "score1")

    eng.train_model(epochs=10, forcenew=False)
    # exit()
    outpath = "Z:\\Group Share\\JGP\\MSV000090488\\short.pkl"
    path = "Z:\\Group Share\\JGP\\MSV000090488\\20220207_UreaExtracted_SW480_C4_RPLC_01.pkl"
    dataset = np.load(path, allow_pickle=True)

    zvals = []
    scores = []
    dummyzs = []
    dummy_scores = []
    correct = 0
    dummycorrect = 0
    for d in dataset:
        centroids = d[0]
        z = d[1]
        score = eng.predictor(centroids, z)
        zvals.append(z)
        scores.append(score)

        if score > 0.5:
            correct += 1

        dummyz = get_dummy_charge_state(z)
        dummy_score = eng.predictor(centroids, dummyz)

        dummyzs.append(dummyz)
        dummy_scores.append(dummy_score)

        if dummy_score < 0.5:
            dummycorrect += 1

    print("Correct Percent:", correct / len(scores))
    print("Dummy Correct Percent:", dummycorrect/ len(dummy_scores), len(dummy_scores))

    # 2D histogram in log scale
    # import lognorm
    from matplotlib.colors import LogNorm
    plt.subplot(221)
    plt.title("True Charge State")
    plt.hist2d(zvals, scores, bins=(np.arange(0.5, 50.5, 1), np.arange(0, 1.1, 0.1)))#, norm=LogNorm())
    plt.ylabel("Score")
    plt.xlabel("Charge")

    plt.subplot(222)
    plt.title("Dummy Charge State")
    plt.hist2d(dummyzs, dummy_scores, bins=(np.arange(0.5, 50.5, 1), np.arange(0, 1.1, 0.1)))#, norm=LogNorm())
    plt.ylabel("Score")
    plt.xlabel("Charge")

    plt.subplot(223)
    # Histograms of scores for both
    plt.hist(scores, bins=np.arange(0, 1.1, 0.1), label="True")
    plt.hist(dummy_scores, bins=np.arange(0, 1.1, 0.1), label="Dummy")
    plt.xlabel("Score")
    plt.ylabel("Count")
    plt.legend()

    plt.subplot(224)
    # Create ROC curve
    from sklearn.metrics import roc_curve, roc_auc_score
    y = np.array([1] * len(scores) + [0] * len(dummy_scores))
    scores = np.array(scores + dummy_scores)
    fpr, tpr, thresholds = roc_curve(y, scores)
    auc = roc_auc_score(y, scores)
    plt.plot(fpr, tpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve, AUC: {auc:.3f}")
    plt.tight_layout()


    plt.show()
