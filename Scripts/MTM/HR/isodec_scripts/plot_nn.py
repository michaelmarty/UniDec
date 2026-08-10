import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import torch
import torch.nn
from unidec.IsoDec.encoding import *

# Set backend to Agg
mpl.use('WxAgg')
os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")

modelpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\modelparams\\phase_model_4.pth"

# modelpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\modelparams\\phase_model_8.pth"

model = torch.load(modelpath, weights_only=True)

w1 = model['linear_relu_stack.0.weight'].cpu().numpy()
b1 = model['linear_relu_stack.0.bias'].cpu().numpy()
w2 = model['linear_relu_stack.2.weight'].cpu().numpy()
b2 = model['linear_relu_stack.2.bias'].cpu().numpy()

print(len(w1.flatten()), len(w2.flatten()), len(b1), len(b2))
print(np.sum([len(w1.flatten()), len(w2.flatten()), len(b1), len(b2)]))
exit()

emat = np.loadtxt("emat.txt")

flat = emat.flatten()

l1 = np.dot(flat, np.transpose(w1)) + b1

l1 = np.maximum(l1, 0)

l2 = np.dot(l1, np.transpose(w2)) + b2

plt.imshow(np.transpose([l2]), aspect=0.25, cmap="viridis", interpolation="none")
plt.show()

exit()

# Plot L2
plt.figure()
plt.imshow(np.transpose([l2]), aspect=0.25, cmap="viridis", interpolation="none")
plt.xticks([])
plt.yticks([])
plt.savefig("l2.png", dpi=300, transparent=True)
plt.savefig("l2.pdf")
plt.show()

exit()
# Plot L1
plt.figure()
plt.imshow(np.transpose([l1]), aspect=0.1, cmap="viridis", interpolation="none")
plt.xticks([])
plt.yticks([])
plt.savefig("l1.png", dpi=300, transparent=True)
plt.savefig("l1.pdf")
plt.show()

exit()



# Plot B1
plt.figure()
plt.imshow(np.transpose([b2]), aspect=0.25, cmap="viridis", interpolation="none")
plt.xticks([])
plt.yticks([])
plt.savefig("b2.png", dpi=300, transparent=True)
plt.savefig("b2.pdf")
plt.show()


plt.figure()
plt.imshow(w2, aspect=1, cmap="viridis", interpolation="none")
plt.xticks([])
plt.yticks([])
plt.savefig("w2.png", dpi=300, transparent=True)
plt.savefig("w2.pdf")
plt.show()
exit()


plt.figure()
plt.imshow(np.transpose([flat]), aspect=0.1, cmap="viridis", interpolation="none")
plt.xticks([])
plt.yticks([])
plt.savefig("inputvector.png", dpi=300, transparent=True)
plt.savefig("inputvector.pdf")
plt.show()
