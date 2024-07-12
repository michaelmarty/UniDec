import numpy as np
import os
import torch
from torch import nn
from torchvision.transforms import ToTensor
from torchvision.models import resnet50, mobilenet_v3_large, MobileNet_V3_Large_Weights, ResNet50_Weights
from torch.optim import lr_scheduler
import time
import platform
from unidec.IsoDec.encoding import encode_isodist

example = np.array([[5.66785531e+02, 1.47770838e+06],
                    [5.67057354e+02, 1.54980838e+06],
                    [5.67507468e+02, 5.21600520e+07],
                    [5.67708173e+02, 8.35557760e+07],
                    [5.67908401e+02, 7.28264240e+07],
                    [5.68060254e+02, 1.87337225e+06],
                    [5.68108674e+02, 4.35435520e+07],
                    [5.68239256e+02, 3.88155375e+06],
                    [5.68309390e+02, 2.05468060e+07],
                    [5.68509951e+02, 7.18109250e+06],
                    [5.68707871e+02, 2.30373500e+06],
                    [5.69150563e+02, 1.57598062e+06],
                    [5.69243121e+02, 1.96390440e+07],
                    [5.69334393e+02, 6.82677120e+07],
                    [5.69425337e+02, 1.22867432e+08],
                    [5.69516492e+02, 1.45702336e+08],
                    [5.69607541e+02, 1.20801936e+08],
                    [5.69698595e+02, 1.06786072e+08],
                    [5.69789906e+02, 6.56232960e+07],
                    [5.69881208e+02, 3.41013880e+07],
                    [5.69972168e+02, 1.70930360e+07],
                    [5.70063432e+02, 9.17621100e+06],
                    [5.70699369e+02, 1.96462650e+06]])


class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(16 * 16 * 3, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 50)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


class IsoDecModel:
    def __init__(self, working_dir="C:\\Data\\IsoNN\\"):
        self.device = None
        self.model = None
        self.loss_fn = None
        self.optimizer = None
        self.scheduler = None
        self.working_dir = working_dir
        self.modelid = None
        self.savepath = None
        self.indexes = None

    def get_model(self, modelid):
        self.modelid = modelid
        savename = "nomodel.pth"
        self.savepath = os.path.join(self.working_dir, savename)

    def setup_model(self, modelid=0):
        self.device = (
            "cuda"
            if torch.cuda.is_available()
            else "mps"
            if torch.backends.mps.is_available()
            else "cpu"
        )
        print(f"Using {self.device} device")

        self.get_model(modelid)
        if os.path.isfile(self.savepath):
            self.load_model()
        self.model = self.model.to(self.device)
        print("Loaded Model:", self.model)

    def setup_training(self):
        if self.model is None:
            self.setup_model()
        self.loss_fn = nn.CrossEntropyLoss()
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)

    def load_model(self):
        if os.path.isfile(self.savepath):
            try:
                self.model.load_state_dict(torch.load(self.savepath))
                print("Model loaded:", self.savepath)
            except Exception as e:
                print("Model failed to load:", self.savepath)
                print(e)
        else:
            print("Model not found:", self.savepath)

    def save_model(self):
        torch.save(self.model.state_dict(), self.savepath)
        print("Model saved:", self.savepath)

    def encode(self, centroids, maxlength=16):
        emat, self.indexes = encode_isodist(centroids, maxlen=maxlength)
        return emat

    def train_model(self, dataloader):
        pass

    def evaluate_model(self, dataloader):
        pass

    def predict(self, centroids):
        pass


class IsoDecClassifier(IsoDecModel):
    def __init__(self, working_dir="C:\\Data\\IsoNN\\"):
        super().__init__(working_dir)

    def get_model(self, modelid):
        self.modelid = modelid
        if modelid == 0:
            self.model = NeuralNetwork()
            savename = "exp_custom_nn_model.pth"
        elif modelid == 1:
            self.model = resnet50()  # weights=ResNet50_Weights.DEFAULT)
            savename = "exp_resnet50_model.pth"
        elif modelid == 2:
            self.model = mobilenet_v3_large()
            savename = "exp_mobilenet_v3_large_model.pth"
        else:
            print("Model ID not recognized")
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def train_model(self, dataloader):
        if self.loss_fn is None or self.optimizer is None:
            self.setup_training()

        size = len(dataloader.dataset)
        self.model.train()
        for batch, (X, y) in enumerate(dataloader):
            y = y.long()
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            pred = self.model(X)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % 100 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader):
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for X, y in dataloader:
                y = y.long()
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X)
                test_loss += self.loss_fn(pred, y).item()
                correct += (pred.argmax(1) == y).type(torch.float).sum().item()
        test_loss /= num_batches
        correct /= size
        print(f"Test Error: \n Accuracy: {(100 * correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

    def predict(self, centroids):
        if self.model is None:
            self.setup_model()
        startime = time.perf_counter()

        test_data = self.encode(centroids)
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)  # test_data.flatten().unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predz = predvec[0].argmax(0)
        print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return int(predz)


class NNmulti(nn.Module):
    def __init__(self, size=16, nfeatures=3):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(size * size * 3, 512),
            nn.ReLU(),
            nn.Linear(512, 128),
            nn.ReLU(),
            nn.Linear(128, size * nfeatures),
            nn.Unflatten(1, (size, nfeatures)),
            nn.Softmax(dim=2),
            nn.Flatten(),
            nn.Linear(size * nfeatures, 128),
            nn.ReLU(),
            nn.Linear(128, 512),
            nn.ReLU(),
            nn.Linear(512, size * nfeatures),
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


class IsoDecSegmenter(IsoDecModel):
    def __init__(self, working_dir=None):
        if working_dir is None:
            if platform.system() == "Linux":
                working_dir = "/xdisk/mtmarty/mtmarty/training/"
            else:
                working_dir = "C:\\Data\\IsoNN\\"

        super().__init__(working_dir)
        self.dims = [2, 32]

    def get_model(self, modelid):
        self.modelid = modelid
        if modelid == 0:
            self.model = NNmulti(size=self.dims[1], nfeatures=self.dims[0])
            savename = "exp_custom_nn_multi_model.pth"
        else:
            print("Model ID not recognized")
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def setup_training(self):
        if self.model is None:
            self.setup_model()
        self.loss_fn = nn.MSELoss()
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)

    def train_model(self, dataloader):
        if self.loss_fn is None or self.optimizer is None:
            self.setup_training()

        size = len(dataloader.dataset)
        self.model.train()
        for batch, (X, y) in enumerate(dataloader):
            # flatten last two dimensions of y
            # y = y.long()
            y = y.float()
            # y = torch.flatten(y, start_dim=1)
            X, y = X.to(self.device), y.to(self.device)
            # Compute prediction error
            pred = self.model(X).reshape(y.shape)

            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % 100 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader):
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for X, y in dataloader:
                y = y.float()
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X).reshape(y.shape)
                test_loss += self.loss_fn(pred, y).item()

                pred = torch.round(pred)

                corr = torch.all(pred == y, dim=1)
                corr = torch.all(corr, dim=1)

                correct += corr.type(torch.float).sum().item()
        test_loss /= num_batches
        correct /= size
        print(f"Test Error: \n Accuracy: {(100 * correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

    def predict(self, centroids):
        if self.model is None:
            self.setup_model()
        startime = time.perf_counter()

        test_data = self.encode(centroids, maxlength=self.dims[1])
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            pred = self.model(x).reshape(self.dims)
            pred = torch.round(pred)

        print("Time:", time.perf_counter() - startime)
        return pred.cpu().data.numpy()


if __name__ == "__main__":
    model = IsoDecClassifier()
    z = model.predict(example)
    print(z)

    pass
