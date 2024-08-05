import numpy as np
import os
import torch
from torch import nn
from torchvision.transforms import ToTensor
# from torchvision.models import resnet50, mobilenet_v3_large, MobileNet_V3_Large_Weights, ResNet50_Weights
from torch.optim import lr_scheduler
import time
import platform
# from torchshape import tensorshape
from unidec.IsoDec.encoding import encode_isodist, encode_phase
import inspect

torch.set_float32_matmul_precision("medium")


def set_debug_apis(state: bool = False):
    torch.autograd.profiler.profile(enabled=state)
    torch.autograd.profiler.emit_nvtx(enabled=state)
    torch.autograd.set_detect_anomaly(mode=state)


set_debug_apis(state=False)

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
# Note: Charge states 5 and 11
'''
# Old model never worked well. Keeping it only for reference.
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
'''

'''
class ChargeClassifierNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification.
    Inputs are 16x16x3 images. Output is len 50 array of probabilities for each charge state 0 to  49.
    """

    def __init__(self, size=16, nfeatures=3, outsize=50):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(size * size * nfeatures, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, outsize)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits
'''


def save_model_to_binary(model, outfile):
    params = model.parameters()
    output = []
    for m in params:
        flat = m.flatten()
        output.append(flat)

    output = torch.cat(output).cpu().detach().numpy()

    output.tofile(outfile)

def print_model(model):
    for name, param in model.named_parameters():
        print(name, param.shape)
        print(param[0], param[1], param[2])


class IsoDecModel:
    """
    Generic IsoDec Model base class. Contains methods for training, evaluation, and prediction.
    Inherited by other models. Not meant to run on it's own.
    """

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
        """
        Get the model based on the model ID.
        :param modelid: Model ID
        :return: None
        """
        self.modelid = modelid
        savename = "nomodel.pth"
        self.savepath = os.path.join(self.working_dir, savename)

    def setup_model(self, modelid=None):
        """
        Setup model and load if savepath exists. Set device.
        :param modelid: Model ID passed to self.get_model()
        :return: None
        """
        if modelid is None:
            modelid = self.modelid
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
            print("Loaded Weights:", self.savepath)
        self.model = self.model.to(self.device)
        print("Loaded Model:", self.model)

    def setup_training(self):
        """"
        Setup loss function, optimizer, and scheduler.
        :return: None
        """
        if self.model is None:
            self.setup_model()
        self.loss_fn = nn.CrossEntropyLoss()
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)

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

    def encode(self, centroids, maxlength=16):
        """
        Encode the centroids into a format for the model.
        :param centroids: Centroid data, m/z and intensity
        :param maxlength: Max length of the encoded data. Will generate 3 x maxlength x maxlength array
        :return: Output array of the encoded data, size 3 x maxlength x maxlength
        """
        emat, self.indexes = encode_isodist(centroids, maxlen=maxlength)
        return emat

    def train_model(self, dataloader):
        """
        Unimplemented. To be implemented by child classes.
        :param dataloader: Data Loader from Training Data
        :return: None
        """
        pass

    def evaluate_model(self, dataloader):
        """
        Unimplemented. To be implemented by child classes.
        :param dataloader: Data Loader from Test Data
        :return: None
        """
        pass

    def predict(self, centroids):
        """
        Unimplemented. To be implemented by child classes.
        :param centroids: Centroid data, m/z and intensity
        :return: Prediction, format may vary depending on model
        """
        pass


'''
class IsoDecClassifier(IsoDecModel):
    """
    IsoDec Classifier Model. Inherits from IsoDecModel. Contains methods for training, evaluation, and prediction.
    Used to classify a set of centroids into a charge state, typically 1-50. Uses a neural network.
    Default net is simple ChargeClassifierNeuralNetwork. Can also use more fancy ones though like ResNet50 or MobileNetV3.
    """

    def __init__(self, working_dir="C:\\Data\\IsoNN\\"):
        super().__init__(working_dir)
        self.modelid = 0
        self.maxlen = 16
        self.nfeatures = 3
        self.maxz = 50

    def get_model(self, modelid):
        """
        Get the model based on the model ID.
        :param modelid: Model ID. 0 is custom, 1 is ResNet50, 2 is MobileNetV3
        :return: None
        """
        self.modelid = modelid
        if modelid == 0:
            self.model = ChargeClassifierNeuralNetwork(size=self.maxlen, nfeatures=self.nfeatures, outsize=self.maxz)
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
        for batch, s in enumerate(dataloader):
            X = s["emat"]
            y = s["z"]
            m = s["mask"]
            if False:
                mmatrix = torch.einsum("abi,abj->abij", m, m)
                ones = torch.ones_like(mmatrix)
                mmatrix = torch.cat([ones, ones, mmatrix], dim=1)
                X = X * mmatrix
            X = X.float()
            y = y.long()
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            pred = self.model(X)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % 1000 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader, maxbad=None):
        baddata = []
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for s in dataloader:
                X = s["emat"]
                y = s["z"]
                y = y.long()
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X)
                test_loss += self.loss_fn(pred, y).item()

                pred_class = pred.argmax(1)
                ccorr = (pred_class == y).type(torch.float)

                correct += ccorr.sum().item()

                if maxbad is not None and len(baddata) < maxbad:
                    for i in range(len(ccorr)):
                        if not ccorr[i]:
                            baddata.append((X[i], s["mask"], y[i], s["mask"], pred_class[i]))

        test_loss /= num_batches
        correct /= size
        print(f"Test Error: \n Accuracy: {(100 * correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

        return baddata

    def predict(self, centroids):
        if self.model is None:
            self.setup_model()
        # startime = time.perf_counter()

        test_data = self.encode(centroids, maxlength=self.maxlen)
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)  # test_data.flatten().unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predz = predvec[0].argmax(0)
        # print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return int(predz), np.ones((1, len(centroids)))


class SegmentationNeuralNet(nn.Module):
    def __init__(self, size=16, nfeatures=3):
        super().__init__()
        self.flatten = nn.Flatten()
        k1 = 3
        k2 = 3
        h1 = 3
        h2 = 16
        self.conv1 = nn.Conv2d(3, h1, k1)
        self.conv2 = nn.Conv2d(h1, h2, k2)
        cout = tensorshape(self.conv2, (tensorshape(self.conv1, (1, 3, size, size))))
        ncout = cout[1] * cout[2] * cout[3]
        self.linear_relu_stack = nn.Sequential(
            self.conv1, nn.ReLU(), self.conv2, nn.ReLU(),
            nn.Flatten(),

            nn.Linear(ncout, 512),
            nn.ReLU(),

            nn.Linear(512, size * nfeatures),
            nn.ReLU(),
        )

    def forward(self, x):
        # x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return torch.clamp_max(logits, 1)


class IsoDecSegmenter(IsoDecModel):
    def __init__(self, working_dir=None):
        if working_dir is None:
            if platform.system() == "Linux":
                working_dir = "/xdisk/mtmarty/mtmarty/training/"
            else:
                working_dir = "C:\\Data\\IsoNN\\"

        super().__init__(working_dir)
        self.dims = [1, 16]
        self.modelid = 0

    def get_model(self, modelid):
        self.modelid = modelid
        if modelid == 0:
            self.model = SegmentationNeuralNet(size=self.dims[1], nfeatures=self.dims[0])
            savename = "exp_custom_nn_multi_model2.pth"
        else:
            print("Model ID not recognized", modelid)
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def setup_training(self):
        if self.model is None:
            self.setup_model(self.modelid)
        self.loss_fn = nn.MSELoss()
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)

    def train_model(self, dataloader):
        if self.loss_fn is None or self.optimizer is None:
            self.setup_training()

        size = len(dataloader.dataset)
        self.model.train()
        for batch, s in enumerate(dataloader):
            X = s["emat"]
            y = s["mask"]
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

            if batch % 1000 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader, maxbad=None):
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct, corr1, corr2, norm2 = 0, 0, 0, 0, 0
        with torch.no_grad():
            for s in dataloader:
                X = s["emat"]
                y = s["mask"]
                y = y.float()
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X).reshape(y.shape)

                test_loss += self.loss_fn(pred, y).item()

                pred = torch.round(pred)
                c1 = pred == y
                c2 = y == 1
                c3 = c1 & c2
                corr1 += torch.sum(c1).item()
                corr2 += torch.sum(c3).item()
                norm2 += torch.sum(c2).item()
                corr = torch.all(c1, dim=1)
                corr = torch.all(corr, dim=1)

                correct += corr.type(torch.float).sum().item()
        test_loss /= num_batches
        correct /= size
        c1 = corr1 / (size * self.dims[0] * self.dims[1])
        c2 = corr2 / norm2
        print(f"Test Error: \n Full Accuracy: {(100 * correct):>0.1f}%, "
              f"Partial Accuracy: {(100 * c1):>0.1f}%, "
              f"Mask Accuracy: {(100 * c2):>0.1f}%, "
              f"Avg loss: {test_loss:>8f} \n")

        return []

    def predict(self, centroids):
        """
        Predict the mask for the input centroids.
        :param centroids: Centroid data, m/z and intensity
        :return: Predicted mask, same length as centroids
        """
        if self.model is None:
            self.setup_model()
        # startime = time.perf_counter()

        test_data = self.encode(centroids, maxlength=self.dims[1])
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            pred = self.model(x).reshape(self.dims)
            pred = torch.round(pred)

        # print("Time:", time.perf_counter() - startime)

        ind = self.indexes[:self.dims[1]]
        mask = pred.cpu().data.numpy().flatten().astype(bool)
        maskedindexes = ind[mask].astype(int)

        return 0, maskedindexes




class CombinedSegClassNeuralNet(nn.Module):
    def __init__(self, size=16, nfeatures=3, outsize=50):
        super().__init__()
        self.flatten = nn.Flatten()
        k1 = 3
        k2 = 3
        h1 = 3
        h2 = 16
        self.conv1 = nn.Conv2d(3, h1, k1)
        self.conv2 = nn.Conv2d(h1, h2, k2)
        cout = tensorshape(self.conv2, (tensorshape(self.conv1, (1, 3, size, size))))
        ncout = cout[1] * cout[2] * cout[3]
        self.seg_stack = nn.Sequential(
            self.conv1, nn.ReLU(), self.conv2, nn.ReLU(),
            nn.Flatten(),

            nn.Linear(ncout, 512),
            nn.ReLU(),

            nn.Linear(512, size * nfeatures),
            nn.ReLU(),
        )

        self.class_stack = nn.Sequential(
            nn.Linear(size * size * 3, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, outsize)
        )

    def forward(self, x):
        # Get the segmentation
        seglogits = torch.clamp_max(self.seg_stack(x), 1)
        # Mask the input for the classifier based on the segmentation output
        mvec = torch.round(seglogits)
        mmatrix = torch.einsum("bi,bj->bij", mvec, mvec)
        cinput = torch.einsum("aij,abij->abij", mmatrix, x)
        # Run the classifier
        classlogits = self.class_stack(self.flatten(cinput))
        return seglogits, classlogits


class CustomLoss(nn.Module):
    def __init__(self):
        super(CustomLoss, self).__init__()
        self.loss_seg = nn.MSELoss()
        self.loss_class = nn.CrossEntropyLoss()

    def forward(self, pred_seg, y, pred_class, z):
        loss_seg = self.loss_seg(pred_seg, y)
        loss_class = self.loss_class(pred_class, z)
        return loss_seg ** 2 + loss_class ** 2


class IsoDecMixedModel(IsoDecModel):
    def __init__(self, working_dir=None):
        if working_dir is None:
            if platform.system() == "Linux":
                working_dir = "/xdisk/mtmarty/mtmarty/training/"
            else:
                working_dir = "C:\\Data\\IsoNN\\"

        super().__init__(working_dir)
        self.dims = [1, 16]
        self.modelid = 1

    def get_model(self, modelid):
        self.modelid = modelid
        if modelid == 0:
            self.model = CombinedSegClassNeuralNet(size=self.dims[1], nfeatures=self.dims[0])
            savename = "combined_model.pth"
        elif modelid == 1:
            self.model = CombinedSegClassNeuralNet(size=self.dims[1], nfeatures=self.dims[0])
            savename = "combined_model_2.pth"
        elif modelid == 2:
            self.model = CombinedSegClassNeuralNet(size=self.dims[1], nfeatures=self.dims[0])
            savename = "combined_model_3.pth"
        else:
            print("Model ID not recognized", modelid)
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def setup_training(self):
        if self.model is None:
            self.setup_model(self.modelid)
        self.loss_fn = CustomLoss()
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        # self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1, weight_decay=0.0001, momentum=0.9)
        # self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.01, weight_decay=0.0001)
        self.scheduler = lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.95)

    def train_model(self, dataloader):
        if self.loss_fn is None or self.optimizer is None:
            self.setup_training()

        size = len(dataloader.dataset)
        self.model.train()
        for batch, s in enumerate(dataloader):
            X = s["emat"]
            y = s["mask"]
            z = s["z"]

            # Load All
            z = z.long()
            y = y.float()
            X, y, z = X.to(self.device), y.to(self.device), z.to(self.device)
            # Predict
            pred_seg, pred_class = self.model(X)
            pred_seg = pred_seg.reshape(y.shape)
            # Compute prediction error
            loss = self.loss_fn(pred_seg, y, pred_class, z)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()
            # Report
            if batch % 1000 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader, maxbad=None):
        if self.model is None:
            self.setup_training()
        baddata = []
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct, corr1, corr2, norm2, classcorr = 0, 0, 0, 0, 0, 0
        with torch.no_grad():
            for s in dataloader:
                X = s["emat"]
                y = s["mask"]
                z = s["z"]
                z = z.long()
                y = y.float()
                X, y, z = X.to(self.device), y.to(self.device), z.to(self.device)
                # Compute prediction error
                pred_seg, pred_class = self.model(X)
                pred_seg = pred_seg.reshape(y.shape)

                test_loss += self.loss_fn(pred_seg, y, pred_class, z).item()

                pred_seg = torch.round(pred_seg)
                c1 = pred_seg == y
                c2 = y == 1
                c3 = c1 & c2
                corr1 += torch.sum(c1).item()
                corr2 += torch.sum(c3).item()
                norm2 += torch.sum(c2).item()
                corr = torch.all(c1, dim=1)
                corr = torch.all(corr, dim=1)
                correct += corr.type(torch.float).sum().item()

                pred_class = pred_class.argmax(1)
                ccorr = (pred_class == z).type(torch.float)
                classcorr += ccorr.sum().item()

                if maxbad is not None and len(baddata) < maxbad:
                    for i in range(len(ccorr)):
                        if not ccorr[i]:
                            baddata.append((X[i], y[i], z[i], pred_seg[i], pred_class[i]))

        test_loss /= num_batches
        correct /= size
        c1 = corr1 / (size * self.dims[0] * self.dims[1])
        c2 = corr2 / norm2

        classcorr /= size

        print(f"Test Error: \n Full Accuracy: {(100 * correct):>0.1f}%, "
              f"Partial Accuracy: {(100 * c1):>0.1f}%, "
              f"Mask Accuracy: {(100 * c2):>0.1f}%, "
              f"Class Accuracy: {(100 * classcorr):>0.1f}%, "
              f"Avg loss: {test_loss:>8f} \n")

        return baddata

    def predict(self, centroids):
        """
        Predict the mask for the input centroids.
        :param centroids: Centroid data, m/z and intensity
        :return: Predicted charge (int), Predicted mask (bool array, same length as centroids)
        """
        if self.model is None:
            self.setup_model()
        startime = time.perf_counter()

        test_data = self.encode(centroids, maxlength=self.dims[1])
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            pred_seg, pred_charge = self.model(x)
            pred_seg = torch.round(pred_seg)
            pred_seg = pred_seg.cpu().data.numpy()

            predz = int(pred_charge[0].argmax(0))

        v = pred_seg[0].astype(bool)
        if len(self.indexes) > len(v):
            v = np.append(v, np.zeros(len(self.indexes) - len(v)))
        else:
            v = v[:len(self.indexes)]
        i = self.indexes[v]
        mask = np.zeros(len(centroids))
        mask[self.indexes[v].astype(int)] = 1
        mask = mask.astype(bool)

        print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return predz, mask
'''


class PhaseNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification.
    Inputs are 50x16 phase images. Output is len 50 array of probabilities for each charge state 0 to  49.
    """

    def __init__(self, size=8, outsize=50):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(size * outsize, size * outsize),
            nn.ReLU(),
            nn.Linear(size * outsize, outsize)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


class FastPhaseNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification.
    Inputs are 50x16 phase images. Output is len 50 array of probabilities for each charge state 0 to  49.
    """

    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(400, 400),
            nn.ReLU(),
            # nn.Linear(400, 400),
            # nn.ReLU(),
            nn.Linear(400, 50)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


class PhaseModel(IsoDecModel):
    def __init__(self, working_dir=None):
        if working_dir is None:
            if platform.system() == "Linux":
                # working_dir = "/xdisk/mtmarty/mtmarty/training/"
                filename = inspect.getframeinfo(inspect.currentframe()).filename
                working_dir = os.path.dirname(os.path.abspath(filename))
            else:
                # working_dir = "C:\\Data\\IsoNN\\"
                filename = inspect.getframeinfo(inspect.currentframe()).filename
                working_dir = os.path.dirname(os.path.abspath(filename))

        super().__init__(working_dir)
        self.dims = [50, 8]
        self.modelid = 1
        #self.get_model(self.modelid)

    def get_model(self, modelid):
        self.modelid = modelid
        if modelid == 0:
            self.model = PhaseNeuralNetwork(size=self.dims[1], outsize=self.dims[0])
            savename = "phase_model_2.pth"
        elif modelid == 1:
            self.model = FastPhaseNeuralNetwork()
            # self.model = torch.compile(self.model, mode="max-autotune")
            savename = "phase_model.pth"
        else:
            print("Model ID not recognized", modelid)
            raise ValueError("Model ID not recognized")

        self.savepath = os.path.join(self.working_dir, savename)

    def train_model(self, dataloader):
        # if self.loss_fn is None or self.optimizer is None:
        self.setup_training()
        set_debug_apis(state=False)
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.train()

        for batch, s in enumerate(dataloader):
            X = s[0]
            y = s[1]
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            pred = self.model(X)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            loss.backward()
            self.optimizer.step()
            self.optimizer.zero_grad()

            if batch % int(num_batches / 5) == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        if self.scheduler is not None:
            self.scheduler.step()

    def evaluate_model(self, dataloader, savebad=False):
        if self.model is None:
            self.setup_training()

        baddata = []
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for j, s in enumerate(dataloader):
                X = s[0]
                y = s[1]
                X, y = X.to(self.device), y.to(self.device)

                pred = self.model(X)
                test_loss += self.loss_fn(pred, y).item()

                pred_class = pred.argmax(1)
                ccorr = (pred_class == y).type(torch.float)

                correct += ccorr.sum().item()

                if savebad:
                    for i in range(len(ccorr)):
                        if not ccorr[i]:
                            baddata.append([j, i, int(y[i]), int(pred_class[i])])

        test_loss /= num_batches
        correct /= size
        print(f"Test Error: \n Accuracy: {(100 * correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

        return baddata

    def predict(self, centroids):
        if self.model is None:
            self.setup_model()
        # startime = time.perf_counter()

        test_data = self.encode(centroids)
        test_data = torch.tensor(test_data, dtype=torch.float32)

        self.model.eval()

        x = test_data.unsqueeze(0)
        with torch.no_grad():
            x = x.to(self.device)
            predvec = self.model(x)
            predz = predvec[0].argmax(0)
        # print("Z:", int(predz), "Time:", time.perf_counter() - startime)
        return int(predz)

    def batch_predict(self, dataloader):
        if self.model is None:
            self.setup_training()
        size = len(dataloader.dataset)
        self.model.eval()
        output = torch.zeros(size, dtype=torch.long, device=self.device)
        with torch.no_grad():
            for batch, x in enumerate(dataloader):
                x = x.to(self.device)
                predvec = self.model(x)
                predz = predvec.argmax(dim=1)
                lx = len(x)
                start = batch * lx
                end = batch * lx + lx
                output[start:end] = predz
        output = output.cpu().numpy()
        return output

    def encode(self, centroids, maxlen=None):
        """
        Encode the centroids into a format for the model.
        :param centroids: Centroid data, m/z and intensity
        :param maxlen: Unused
        :return: Output array of the encoded data, size 3 x maxlength x maxlength
        """
        emat = encode_phase(centroids, maxz=self.dims[0], phaseres=self.dims[1])
        return emat


if __name__ == "__main__":
    model = PhaseModel()
    model.save_model()
    z = model.predict(example)
    print(z)

    pass
