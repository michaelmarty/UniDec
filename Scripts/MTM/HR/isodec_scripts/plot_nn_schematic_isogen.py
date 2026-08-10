import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import unidec.IsoDec.isogenwrapper as ig

os.chdir("C:\\Users\\mm96978\\The University of Texas at Austin\\MartyLab - General\\Papers\\IsoGen\\Figures")
# os.chdir("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\")
# Create plot to show two linear layers connected, the first with 200 nodes and the second with 50

# Set backend to Agg
mpl.use('WxAgg')

if True:
    alpha = 0.2

    n1 = 5
    n2 = 32
    s1 = (n2-1)/(n1-1)

    for n in range(0, n1):
        for m in range(0, n2):
            plt.plot([0, 1], [n * s1, m], color="royalblue", alpha=alpha, marker="o", markersize=1)

    for n in range(0, n2):
        for m in range(0, n2):
            plt.plot([1, 2], [n, m], color="b", alpha=alpha, marker="o", markersize=1)

    for n in range(0, n2):
        for m in range(0, n2):
            plt.plot([2, 3], [n, m], color="darkblue", alpha=alpha, marker="o", markersize=1)

    # Turn off all axes
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("nn_schematic_5_16.png", dpi=300, transparent=True)
    plt.savefig("nn_schematic_5_16.pdf", dpi=300, transparent=True)
    plt.show()

if False:
    mass = 2000
    isodist = ig.nn_gen_isodist(mass, isolen=32)
    color1 = "indigo"

    # Set figure to 0.5 inches wide and 1.25 inches tall
    plt.figure(figsize = (0.5, 1.25))
    # Plot isodist as a stick plot for the first 5 isotopes
    for i in range(7):
        plt.vlines(i, 0, isodist[i], color=color1, alpha=0.8, linewidth=2)
        plt.plot(i, isodist[i], marker="o", color=color1, markersize=3)

    # Turn off top and right axes. set range from 0 to 1 on y-axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylim(0, 1.05)
    plt.xlim(-0.5, 7.5)
    # turn off all axis labels and ticks
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)
    plt.xticks([])
    plt.yticks([])

    # Set line thickness on axes to 2
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)


    # Set aspect ratio to be tall and skinny
    plt.gca().set_aspect(20)
    plt.tight_layout()
    plt.savefig("nn_schematic_isodist.png", dpi=300, transparent=True)
    plt.savefig("nn_schematic_isodist.pdf", dpi=300, transparent=True)
    plt.show()
