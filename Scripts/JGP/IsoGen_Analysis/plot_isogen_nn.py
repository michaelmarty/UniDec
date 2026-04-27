import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import torch
import torch.nn
from unidec.IsoDec.IsoGen.isogen_tools import mass_to_vector
from matplotlib import rcParams


def softsign(x):
    return x / (1 + np.abs(x))

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Set backend to Agg
mpl.use('WxAgg')
os.chdir(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures")

rcParams['axes.linewidth'] = 1
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

modelpath = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\OriginalModels\isogenmass_model_32.pth"


model = torch.load(modelpath, weights_only=True)

w1 = model['linear_relu_stack.0.weight'].cpu().numpy()
b1 = model['linear_relu_stack.0.bias'].cpu().numpy()
w2 = model['linear_relu_stack.2.weight'].cpu().numpy()
b2 = model['linear_relu_stack.2.bias'].cpu().numpy()
w3 = model['linear_relu_stack.4.weight'].cpu().numpy()
b3 = model['linear_relu_stack.4.bias'].cpu().numpy()

mass = 9268.56 #Da
vec = mass_to_vector(mass)

flat = vec.flatten()

h1 = np.dot(flat, np.transpose(w1)) + b1

#Softsign transformation
h1 = np.array([softsign(x) for x in h1])

h2 = np.dot(h1, np.transpose(w2)) + b2

h2 = np.array([softsign(x) for x in h2])

out = np.dot(h2, np.transpose(w3)) + b3

out = np.array([sigmoid(x) for x in out])


# #plot the input vector
# print("Plotting input vector")
# plt.figure()
# plt.imshow(np.transpose([flat]), aspect=3, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("inputvector.png", dpi=300, transparent=True)
# plt.savefig("inputvector.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot w1
# print("Plotting w1")
# plt.figure()
# plt.imshow(np.transpose([w1]), aspect=20, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("w1.png", dpi=300, transparent=True)
# plt.savefig("w1.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot b1
# print("Plotting b1")
# plt.figure()
# plt.imshow(np.transpose([b1]), aspect=0.25, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("b1.png", dpi=300, transparent=True)
# plt.savefig("b1.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot h1
# print("Plotting h1")
# plt.figure()
# plt.imshow(np.transpose([h1]), aspect=0.25, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("h1.png", dpi=300, transparent=True)
# plt.savefig("h1.svg", dpi=300, transparent=True)
# plt.show()


# # Plot w2
# print("Plotting w2")
# plt.figure()
# plt.imshow(np.transpose([w2]), aspect=1, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("w2.png", dpi=300, transparent=True)
# plt.savefig("w2.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot b2
# print("Plotting b2")
# plt.figure()
# plt.imshow(np.transpose([b2]), aspect=0.25, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("b2.png", dpi=300, transparent=True)
# plt.savefig("b2.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot h2
# print("Plotting h2")
# plt.figure()
# plt.imshow(np.transpose([h2]), aspect=0.25, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("h2.png", dpi=300, transparent=True)
# plt.savefig("h2.svg", dpi=300, transparent=True)
# plt.show()


# # Plot w3
# print("Plotting w3")
# plt.figure()
# plt.imshow(np.transpose([w3]), aspect=1, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("w3.png", dpi=300, transparent=True)
# plt.savefig("w3.svg", dpi=300, transparent=True)
# plt.show()
#
#
# # Plot b3
# print("Plotting b3")
# plt.figure()
# plt.imshow(np.transpose([b3]), aspect=0.25, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("b3.png", dpi=300, transparent=True)
# plt.savefig("b3.svg", dpi=300, transparent=True)
# plt.show()


# Plot out
# print("Plotting out")
# plt.figure()
# plt.imshow(np.transpose([out]), aspect=1, cmap="viridis", interpolation="none")
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.xticks([])
# plt.yticks([])
# plt.savefig("out.png", dpi=300, transparent=True)
# plt.savefig("out.svg", dpi=300, transparent=True)
# plt.show()

plt.figure()
plt.plot(out, color='k')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xticks(np.arange(32))
plt.tick_params(axis="both", length=5, width=1, labelbottom=False)
plt.yticks([0,1], fontsize=16)
plt.xlim(0,31)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig("out_plot.png", dpi=300, transparent=True)
plt.savefig("out_plot.pdf", dpi=300, transparent=True)
plt.show()

