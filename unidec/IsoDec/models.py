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
from unidec.IsoDec.encoding import  encode_phase
import inspect

# Speed up matmul precision for faster training
torch.set_float32_matmul_precision("medium")


# Turn on debug APIs to speed up training
def set_debug_apis(state: bool = False):
    torch.autograd.profiler.profile(enabled=state)
    torch.autograd.profiler.emit_nvtx(enabled=state)
    torch.autograd.set_detect_anomaly(mode=state)


set_debug_apis(state=False)

# Example data
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

# Function to save the model to a binary format that can be read into C
def save_model_to_binary(model, outfile):
    params = model.parameters()
    output = []
    for m in params:
        flat = m.flatten()
        output.append(flat)

    output = torch.cat(output).cpu().detach().numpy()

    output.tofile(outfile)


# Function to print the model parameters
def print_model(model):
    for name, param in model.named_parameters():
        print(name, param.shape)
        print(param[0], param[1], param[2])


class PhaseModel:
    """
    General model class for charge state prediction base on a phase encoding.

    Includes functions to train, evaluate, and predict charge states.
    """

    def __init__(self, working_dir=None):
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
        self.dims = [50, 8]
        self.class_weights = None

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
            self.model = Fast4PhaseNeuralNetwork()
            self.dims = [50, 4]
            savename = "phase_model_4.pth"
        elif modelid == 1:
            self.model = Fast8PhaseNeuralNetwork()
            self.dims = [50, 8]
            # self.model = torch.compile(self.model, mode="max-autotune")
            savename = "phase_model_8.pth"
        elif modelid == 2:
            self.model = PhaseNeuralNetwork(size=self.dims[1], outsize=self.dims[0])
            savename = "phase_model_3.pth"
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
            else "mps"
            if torch.backends.mps.is_available()
            else "cpu"
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

    def get_class_weights(self, dataloader):
        """
        Calculate class weights for the data in a DataLoader object.
        :param dataloader: Training DataLoader object
        :return: None
        """
        # Get all the charge assignments
        zvals = torch.cat([x[1] for x in dataloader])
        # Histogram the values
        weights = torch.histogram(zvals.float(), bins=torch.arange(0, 51, dtype=torch.float32))[0]

        # Invert the histogram and normalize
        b1 = weights > 0
        weights[b1] = 1.0 / (weights[b1])
        weights = weights / torch.max(weights)

        # Set the class weights
        self.class_weights = weights
        # self.class_weights.to(self.device)
        print("Class Weights:", self.class_weights, len(self.class_weights))

    def setup_training(self, lossfn="crossentropy", forcenew=False):
        """"
        Setup loss function, optimizer, and scheduler.
        :param lossfn: Loss function to use. Options are "crossentropy", "weightedcrossentropy", and "focal".
        :return: None
        """
        if self.model is None:
            self.setup_model(forcenew=forcenew)
        if lossfn == "crossentropy":
            self.loss_fn = nn.CrossEntropyLoss()
        elif lossfn == "weightedcrossentropy":
            if self.class_weights is None:
                raise ValueError("Class weights not set.")
            self.loss_fn = nn.CrossEntropyLoss(weight=self.class_weights.to(self.device))
        elif lossfn == "focal":
            if self.class_weights is None:
                raise ValueError("Class weights not set.")
            self.loss_fn = FocalLoss(alpha=self.class_weights.to(self.device))
        else:
            raise ValueError("Loss function not recognized.", lossfn)
        # self.loss_fn = nn.CrossEntropyLoss()
        # self.loss_fn = FocalLoss(alpha=self.class_weights.to(self.device))
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

    def train_model(self, dataloader, lossfn="crossentropy", forcenew=False):
        """
        Train the model on a DataLoader object.
        :param dataloader: Training DataLoader object
        :param lossfn: Loss function to use. Options are "crossentropy", "weightedcrossentropy", and "focal".
        :return: None
        """
        # if self.loss_fn is None or self.optimizer is None:
        # plot_zdist(self)
        self.setup_training(lossfn, forcenew=forcenew)
        set_debug_apis(state=False)
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.train()

        for batch, s in enumerate(dataloader):
            x = s[0]
            y = s[1]
            x, y = x.to(self.device), y.to(self.device)

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

    def evaluate_model(self, dataloader, savebad=False):
        """
        Evaluate the model on a test set.
        :param dataloader: Test DataLoader object
        :param savebad: Whether to collect incorrect predictions
        :return: List of incorrect predictions if savebad is True
        """
        if self.model is None:
            self.setup_training()

        baddata = []
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for j, s in enumerate(dataloader):
                x = s[0]
                y = s[1]
                x, y = x.to(self.device), y.to(self.device)

                pred = self.model(x)
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
        """
        Predict charge state for a single set of centroids.
        :param centroids: Centroid data, m/z and intensity
        :return: Predicted charge state, integer
        """
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

    def batch_predict(self, dataloader, zscore_thresh=0.95):
        """
        Predict charge states for a batch of data.
        :param dataloader: DataLoader object with the data to predict
        :return: Output charge state predictions
        """
        if self.model is None:
            self.setup_training()
        size = len(dataloader.dataset)
        self.model.eval()
        output = torch.zeros((size, 2), dtype=torch.long, device=self.device)
        with torch.no_grad():
            for batch, x in enumerate(dataloader):
                x = x.to(self.device)
                lx = len(x)
                start = batch * lx
                end = batch * lx + lx
                predvec = self.model(x)
                #print(predvec)
                predzs, predz_inds = torch.topk(predvec, k=2, dim=1)
                output[start:end, 0] = predz_inds[:, 0]
                second_score_within = ((predzs[:, 1] / predzs[:, 0]) > zscore_thresh).float()
                output[start:end, 1] = predz_inds[:, 1] * second_score_within
        output = output.cpu().numpy()
        return output

    def encode(self, centroids):
        """
        Encode the centroids into a format for the model.
        :param centroids: Centroid data, m/z and intensity
        :return: Output array of the encoded data, size dims[0] x dims[1]
        """
        emat = encode_phase(centroids, maxz=self.dims[0], phaseres=self.dims[1])
        return emat


class PhaseNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification. Generalized dimensions.
    Inputs are nxm phase images. Output is len outsize array of probabilities for each charge state 0 to  outsize-1.
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


class Fast8PhaseNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification.
    Inputs are 50x8 phase images. Output is len 50 array of probabilities for each charge state 0 to  49.
    """

    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(400, 400),
            nn.ReLU(),
            nn.Linear(400, 50)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

    '''
    def partial(self, x):
        x = self.flatten(x)

        weights1 = self.linear_relu_stack[0].weight
        h1 = x @ weights1.t()
        b1 = self.linear_relu_stack[0].bias
        h1 = h1 + b1
        h1 = torch.relu(h1)
        h2 = h1 @ self.linear_relu_stack[2].weight.t()
        b2 = self.linear_relu_stack[2].bias
        h2 = h2 + b2
        print(h1[:5])
        print(h2[:5])

        logits = self.linear_relu_stack[0](x)
        logits = self.linear_relu_stack[1](logits)
        logits = self.linear_relu_stack[2](logits)
        return logits'''


class Fast4PhaseNeuralNetwork(nn.Module):
    """
    Very simple neural net for classification.
    Inputs are 50x4 phase images. Output is len 50 array of probabilities for each charge state 0 to  49.
    """

    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(200, 200),
            nn.ReLU(),
            nn.Linear(200, 50)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


class FocalLoss(nn.Module):
    """
    Class for focal loss function.
    """

    def __init__(self, alpha, gamma=2):
        """
        Initialize the focal loss function.
        :param alpha: Class weights. Should be a tensor.
        :param gamma: Power parameter. Default is 2.
        """
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma

    def forward(self, inputs, targets):
        ce_loss = nn.functional.cross_entropy(inputs, targets, reduction='none')
        pt = torch.exp(-ce_loss)
        loss = (self.alpha[targets] * (1 - pt) ** self.gamma * ce_loss).mean()
        return loss


if __name__ == "__main__":
    model = PhaseModel()
    model.save_model()
    z = model.predict(example)
    print(z)

    pass
