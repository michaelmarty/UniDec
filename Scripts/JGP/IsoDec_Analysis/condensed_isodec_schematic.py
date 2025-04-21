import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from unidec.modules.isotopetools import *
from unidec.IsoDec.plots import *
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import *
from matplotlib import rcParams
from unidec.IsoDec.models import *

def get_full_emat(mock_data, phaseres=4):
    emat = encode_phase(mock_data, maxz=50, phaseres=phaseres)
    return emat

def generate_mock_data(masses, charges, intensities):
    distributions = []
    min_mz = 100000
    max_mz = 0

    for i in range(len(masses)):
        distribution = fast_calc_averagine_isotope_dist(mass=masses[i], charge=charges[i])
        #Scale the distribution by the intensity
        distribution[:,1] = distribution[:,1]*intensities[i]
        distributions.append(distribution)

        dist_min_mz = distribution[0][0]
        dist_max_mz = distribution[-1][0]

        if dist_min_mz < min_mz:
            min_mz = dist_min_mz
        if dist_max_mz > max_mz:
            max_mz = dist_max_mz


    left_padding = 0.2
    right_padding = 0.2
    resolution = 100
    spacing = 1 / resolution


    mock_data = []


    #Add the distribution to the mock data
    for j in range(len(distributions)):
        for i in range(len(distributions[j])):
            mock_data.append(distributions[j][i])

    mock_data = np.array(mock_data)
    mock_data[:,1] = mock_data[:,1]*100

    #Sort the mock data by mz
    mock_data = mock_data[mock_data[:,0].argsort()]

    return mock_data

def get_phase_bin(mz, z, phaseres):
    rescale = mz / mass_diff_c
    phase = (rescale * (z)) % 1
    phase_bin = int(np.floor(phase * phaseres))
    return phase_bin

def plot_wave_under_data(mock_data, query_z, phaseres=4):
    xaxis = np.arange(np.amin(mock_data[:, 0]), np.amax(mock_data[:, 0]), 0.0001)
    yaxis = (np.sin(xaxis * (2*np.pi) * (query_z / 1.0033)) / 2 - 0.5) * 100

    phase_divided_wave = [[] for i in range(phaseres)]
    previous_phase_bin = None
    current_points = []
    for i in range(len(yaxis)):
        phase_bin = get_phase_bin(xaxis[i], query_z, phaseres)
        if phase_bin == previous_phase_bin:
            current_points.append([xaxis[i], yaxis[i]])
        else:
            if len(current_points) > 0:
                phase_divided_wave[previous_phase_bin].append(np.array(current_points))
            current_points = [[xaxis[i], yaxis[i]]]
            previous_phase_bin = phase_bin
    phase_divided_data = [[] for i in range(phaseres)]
    for i in range(len(mock_data)):
        phase_bin = get_phase_bin(mock_data[i][0], query_z, phaseres)
        phase_divided_data[phase_bin].append(mock_data[i])

    colors = ["r", "g", "b", "orange", "c", "m", "y", "k"]

    plt.figure(figsize=(3,5))
    for i in range(len(phase_divided_wave)):
        for j in range(len(phase_divided_wave[i])):
            plt.plot(phase_divided_wave[i][j][:, 0], phase_divided_wave[i][j][:, 1], color=colors[i], linewidth=3)

    for i in range(len(phase_divided_data)):
        if len(phase_divided_data[i]) > 0:
            cplot(np.array(phase_divided_data[i]), color=colors[i])


    mz_range = [np.amin(mock_data[:, 0]), np.amax(mock_data[:, 0])]
    #Plot x ticks on every whole number within the mz_range
    xticks = np.arange(math.ceil(mz_range[0]), math.ceil(mz_range[1]), 1)

    yticks = [0, 50, 100]
    yticklabels = ["0", "%", "100"]

    plt.ylim(-100,100)
    plt.yticks(yticks, yticklabels,fontsize=ytick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.xlim(mz_range[0], mz_range[1])
    plt.xticks(xticks, fontsize=xtick_size)
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\condensed_id_wave_under_dat.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\condensed_id_wave_under_dat.pdf", dpi=300, transparent=True)
    plt.show()

    return

def plot_emat_highlight_z(emat, z, phaseres=4):
    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'

    vmin = np.amin(emat)
    vmax = np.amax(emat)
    #First make an isolated plot of only this z
    plt.figure(figsize=(3,4))
    plt.imshow(np.transpose([emat[z-1]]), aspect=5, cmap="viridis", interpolation="none", vmin=vmin, vmax=vmax)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.tick_params(axis='both', which='both', length=0)
    plt.tight_layout()
    if False:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_z3_encoding.svg", dpi=300, transparent=True, format="svg")
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_z3_encoding.png", dpi=300, transparent=True)
    plt.show()

    # Need to make a plot that highlights the portion of the encoded matrix that corresponds to the given z
    min_index = z - 1

    if True:
        plt.figure(figsize=(3, 4))
        plt.imshow(emat, aspect=0.1, cmap="viridis", interpolation="none", vmin=vmin, vmax=vmax)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        # Create a rectangle that highlights the portion of the matrix that corresponds to the given z
        rect = patches.Rectangle((-0.5, min_index - 0.5), phaseres, 1, linewidth=1.5, edgecolor='r', facecolor='none')
        plt.gca().add_patch(rect)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel("Charge", fontsize=ylabel_size, labelpad=-20)
        # axs[1].tick_params(axis='both', which='both', length=0)
        plt.xlim(-1, phaseres)
        plt.ylim(49.5, -1)
        plt.tight_layout()
        if True:
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_highlight_z3.png", dpi=300, transparent=True)
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_highlight_z3.svg", dpi=300, transparent=True,
                        format="svg")
        plt.show()

def predict_plot_predvec(mock_data, phaseres=4):
    phasemodel = PhaseModel()
    if phaseres == 4:
        phasemodel.modelid = 0
    elif phaseres == 8:
        phasemodel.modelid = 1
    predvec = phasemodel.predict_returnvec(mock_data)
    max_index = np.argmax(predvec)
    print("Predicted charge:", max_index)
    #Create a box around the max index
    rect = patches.Rectangle((-0.5, max_index-0.5), 1, 1, linewidth=1.5, edgecolor='r', facecolor='none')
    plt.figure(figsize=(3, 4))
    plt.imshow(np.transpose([predvec]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().add_patch(rect)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_predvec3.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\id_condensed_predvec3.svg", dpi=300, transparent=True, format="svg")
    plt.show()

# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['font.family'] = 'Arial'

xtick_size = 18
ytick_size = 18
xlabel_size = 24
ylabel_size = 14

masses = [2005, 2998, 3989]
charges = [2, 3, 4]
intensities = [0.55,0.85,0.25]

mock_data = generate_mock_data(masses, charges, intensities)

#Find the maxiumum position of the mock data
maxindex = np.argmax(mock_data[:,1])
maxmz = mock_data[maxindex][0]

padding = 0.25



#Add the rectangles to block out the regions that are not selected to be encoded
min_mz = mock_data[0][0] - padding
max_mz = mock_data[-1][0] + padding

min_encoded_mz = maxmz - 1.05
max_encoded_mz = maxmz + 2.05

if False:
    plt.figure(figsize=(3,5))
    for i in range(len(mock_data)):
        plt.plot([mock_data[i][0], mock_data[i][0]], [0, mock_data[i][1]], color='black')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    #Make rectangles for the non-encoded regions
    rect = patches.Rectangle((min_mz,0), width=min_encoded_mz - min_mz,  height=100, linewidth=1, edgecolor='none', facecolor='b', alpha=0.25)
    plt.gca().add_patch(rect)

    rect2 = patches.Rectangle((max_encoded_mz,0), width=max_mz - max_encoded_mz,  height=100, linewidth=1, edgecolor='none', facecolor='b', alpha=0.25)
    plt.gca().add_patch(rect2)

    plt.plot(maxmz, 90, 'b*', markersize=20)

    plt.xlim([min_mz, max_mz])
    plt.ylim([0, 100])
    plt.yticks(np.arange(0,101,50), labels=["", "%", "100"], fontsize=ytick_size)
    plt.xticks(np.arange(998,1007,2), fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\condensed_id_selection.png", dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\condensed_id_selection.pdf", dpi=300)
    plt.show()

if True:
    #Isolate the mock_data within the encoding region
    encoded_mock_data = mock_data[(mock_data[:,0] > min_encoded_mz) & (mock_data[:,0] < max_encoded_mz)]
    plot_wave_under_data(encoded_mock_data, 3, phaseres=4)
    # emat = get_full_emat(encoded_mock_data, phaseres=4)
    # plot_emat_highlight_z(emat, 3, phaseres=4)
    # predict_plot_predvec(encoded_mock_data, phaseres=4)

