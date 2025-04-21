import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from Scripts.JGP.IsoDec_Analysis import encoding_figure
import numpy as np
from matplotlib import rcParams

# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

os.chdir("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\")

#Need to generate mock data of a single isotope distribution
#Then need to generate the encoded matrix that would be fed into the neural network.
#Then need to generate the neural net schematic with lines connecting each node to the next layer
#Then need to generate the final output prediction vector.

mock_data = encoding_figure.generate_mock_data(2000, 5)
# Find the mz at the maximum intensity in the mock data
apex_mz = mock_data[np.argmax([x[1] for x in mock_data])][0]
min_mz = apex_mz - 1.05
max_mz = apex_mz + 2.05

encoding_figure.plot_wave_under_data(mock_data, 5, phaseres=4)
if True:
    plt.figure(figsize=(1,10))
    for i in range(len(mock_data)):
        plt.plot([mock_data[i][0], mock_data[i][0]], [0, mock_data[i][1]], color="k")
        plt.plot([mock_data[i][0], mock_data[i][0]], [0, -mock_data[i][1]], color="b")
    plt.ylim(-100,100)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    #plt.yticks(np.arange(0, 101, 50), labels=[0, "%", 100], fontsize=24)
    plt.yticks([], labels=[])
    plt.xlabel("m/z", fontsize=32)
    plt.hlines(y=0, xmin=400, xmax=405, color="k")
    plt.xticks(np.arange(0, 2000, 1), fontsize=24)
    plt.xlim(400.65, 402.5)
    plt.gca().yaxis.set_ticks_position('left')
    plt.tight_layout()
    if True:
        plt.savefig("z5_mock_data_matched.png", dpi=300, transparent=True)
        plt.savefig("z5_mock_data_matched.pdf", dpi=300, transparent=True)
    plt.show()
    exit()

emat = encoding_figure.get_full_emat(mock_data, 4)
emat = emat.flatten()

if False:
    plt.figure(figsize=(3, 4))
    plt.imshow(np.transpose([emat]), cmap="viridis", aspect="0.025")
    plt.xticks([])
    plt.yticks([])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.savefig("z5_encoded_matrix.png", dpi=300, transparent=True)
    plt.savefig("z5_encoded_matrix.svg", dpi=300, transparent=True)
    plt.show()

encoding_figure.predict_plot_predvec(mock_data, 4)
