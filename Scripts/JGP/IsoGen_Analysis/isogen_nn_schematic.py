import matplotlib.pyplot as plt
import os
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')


os.chdir(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures")
# Create plot to show two linear layers connected, the first with 200 nodes and the second with 50

# for n in range(0, 200):
#     for m in range(0, 200):
#         plt.plot([0, 1], [n, m], color="k", alpha=0.01, marker="o", markersize=1)
#
# for n in range(0, 200):
#     for m in range(0, 50):
#         plt.plot([1, 2], [n, m * 4 + 1], color="k", alpha=0.01, marker="o", markersize=1)

#Let's do this for the protein mass-based model (isolen 32)
#input length: 5
#first hidden layer: 32
#second hidden layer: 32
#output: 32

input_len = 5
outlen = 32

#Find the midpoints:
l1_mid = (input_len - 1) / 2
l2_mid = (outlen - 1) / 2

diff = l2_mid - l1_mid

for n in range(0, input_len):
    for m in range(0, outlen):
        plt.plot([0, 1], [n + diff, m], color="k", alpha=0.1, marker="o", markersize=1)

for n in range(0, outlen):
    for m in range(0, outlen):
        plt.plot([1, 2], [n, m], color="k", alpha=0.1, marker="o", markersize=1)

for n in range(0, outlen):
    for m in range(0, outlen):
        plt.plot([2, 3], [n, m], color="k", alpha=0.1, marker="o", markersize=1)


# Turn off all axes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().xaxis.set_visible(False)
plt.gca().yaxis.set_visible(False)
plt.xticks([])
plt.yticks([])
plt.savefig("isogen_nn_schematic_5in_32out.png", dpi=300, transparent=True)
plt.savefig("isogen_nn_schematic_5in_32out.pdf", dpi=300, transparent=True)
plt.show()