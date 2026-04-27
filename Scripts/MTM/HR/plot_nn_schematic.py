import matplotlib.pyplot as plt
import os
import matplotlib as mpl



# os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")
#os.chdir("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\")
# Create plot to show two linear layers connected, the first with 200 nodes and the second with 50

# Set backend to Agg
mpl.use('WxAgg')



for n in range(0, 200):
    for m in range(0, 200):
        plt.plot([0, 1], [n, m], color="k", alpha=0.01, marker="o", markersize=1)

for n in range(0, 200):
    for m in range(0, 50):
        plt.plot([1, 2], [n, m * 4 + 1], color="k", alpha=0.01, marker="o", markersize=1)

# Turn off all axes
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().xaxis.set_visible(False)
plt.gca().yaxis.set_visible(False)
plt.xticks([])
plt.yticks([])
# plt.savefig("nn_schematic_200.png", dpi=300, transparent=True)
# plt.savefig("nn_schematic_200.pdf", dpi=300, transparent=True)
plt.show()
