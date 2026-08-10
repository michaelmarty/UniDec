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
from unidec.IsoDec.engine import IsoDecEngine


def get_phase_bin(mz, z, phaseres):
    rescale = mz / mass_diff_c
    phase = (rescale * (z)) % 1
    phase_bin = int(np.floor(phase * phaseres))
    return phase_bin

def generate_mock_data(mass, charge):
    distribution = fast_calc_averagine_isotope_dist(mass=mass, charge=charge)
    min_mz = distribution[0][0]
    max_mz = distribution[-1][0]

    left_padding = 0.2
    right_padding = 0.2
    resolution = 100
    spacing = 1 / resolution
    mz_range = [min_mz-left_padding, max_mz+right_padding]

    # noise_mzs = np.random.uniform(mz_range[0], mz_range[1], size=(int((mz_range[1]-mz_range[0])/spacing), 1))
    # noise_intensities = np.random.normal(0, 0.005, size=(noise_mzs.shape[0], 1))
    # noise_intensities = np.abs(noise_intensities)


    mock_data = []

    #Add the distribution to the mock data
    for i in range(len(distribution)):
        mock_data.append(distribution[i])

    # #Add the noise to the mock data
    # for i in range(len(noise_mzs)):
    #     mock_data.append([noise_mzs[i][0], noise_intensities[i][0]])

    mock_data = np.array(mock_data)
    mock_data[:,1] = mock_data[:,1]*100

    return mock_data

def plot_wave_under_data(mock_data, query_z, phaseres=4):
    xtick_size = 24
    ytick_size = 24
    xlabel_size = 32
    ylabel_size = 24

    xaxis = np.arange(np.min(mock_data[:,0])-0.5, np.max(mock_data[:,0]) + 0.5,0.0001)
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

    # colors = ["r", "g", "b", "orange", "c", "m", "y", "k"]
    colors = ["blue", "orange", "green", "purple"]
    plt.figure(figsize=(1,10))
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
    plt.yticks([])
    plt.xlabel("m/z", style='italic', fontsize=xlabel_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, 0, 2000, color='black', linewidth=1)
    plt.xticks(np.arange(0,2000,1), fontsize=xtick_size)
    plt.xlim(854.8, 856.75)
    plt.tight_layout()
    if False:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\wave_under_data7.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\wave_under_data7.pdf", dpi=300, transparent=True)
    plt.show()

    return

def mock_encode_data(mock_data, z, phaseres=4):
    single_z_encoding = np.zeros(phaseres)
    for i in range(len(mock_data)):
        phase_bin = get_phase_bin(mock_data[i][0], z, phaseres)
        single_z_encoding[phase_bin] += mock_data[i][1]
    return single_z_encoding

def get_full_emat(mock_data, phaseres=4):
    emat = encode_phase(mock_data, maxz=50, phaseres=phaseres)
    return emat

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
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\z1_encoding.svg", dpi=300, transparent=True, format="svg")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\z1_encoding.png", dpi=300, transparent=True)
    plt.show()


    #Need to make a plot that highlights the portion of the encoded matrix that corresponds to the given z
    min_index = z-1

    if True:
        plt.figure(figsize=(3,4))
        plt.imshow(emat, aspect=0.1, cmap="viridis", interpolation="none", vmin=vmin, vmax=vmax)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        #Create a rectangle that highlights the portion of the matrix that corresponds to the given z
        rect = patches.Rectangle((-0.5,min_index-0.5), phaseres, 1, linewidth=1.5, edgecolor='r', facecolor='none')
        plt.gca().add_patch(rect)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel("Charge", fontsize=ylabel_size, labelpad=-20)
        # axs[1].tick_params(axis='both', which='both', length=0)
        plt.xlim(-1, phaseres)
        plt.ylim(49.5,-1)
        plt.tight_layout()
        if False:
            plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\highlight_z7_encoding_new.png", dpi=300, transparent=True)
            plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\highlight_z7_encoding_new.svg", dpi=300, transparent=True, format="svg")
        plt.show()

    flat = emat.flatten()



    plt.figure(figsize=(1,4))
    plt.imshow(np.transpose([flat]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if False:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\flat_emat.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\flat_emat.pdf", dpi=300, transparent=True)
    plt.show()
    return

def predict_plot_predvec(mock_data, phaseres=4):
    phasemodel = PhaseModel()
    if phaseres == 4:
        phasemodel.modelid = 0
    elif phaseres == 8:
        phasemodel.modelid = 1
    predvec = phasemodel.predict_returnvec(mock_data)
    max_index = np.argmax(predvec)
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
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\predvec_z5.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\predvec_z5.svg", dpi=300, transparent=True, format="svg")
    plt.show()

def plot_nn_layers(emat):
    modelpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\modelparams\\phase_model_4.pth"

    # modelpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\modelparams\\phase_model_8.pth"

    model = torch.load(modelpath, weights_only=True)

    w1 = model['linear_relu_stack.0.weight'].cpu().numpy()
    b1 = model['linear_relu_stack.0.bias'].cpu().numpy()
    w2 = model['linear_relu_stack.2.weight'].cpu().numpy()
    b2 = model['linear_relu_stack.2.bias'].cpu().numpy()


    flat = emat.flatten()

    l1 = np.dot(flat, np.transpose(w1)) + b1

    l1 = np.maximum(l1, 0)

    l2 = np.dot(l1, np.transpose(w2)) + b2

    #Plot w1
    plt.figure(figsize=(4,4))
    plt.imshow(w1, aspect=1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\w1_jgp.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\w1_jgp.pdf", dpi=300, transparent=True)
    plt.show()

    #Plot w2
    plt.figure(figsize=(4,4))
    plt.imshow(np.transpose(w2), aspect=1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\w2_jgp.png", dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\w2_jgp.pdf", dpi=300, transparent=True)
    plt.show()

    #Plot b1
    plt.figure(figsize=(1, 4))
    plt.imshow(np.transpose([b1]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\b1_jgp.png",
                    dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\b1_jgp.pdf",
                    dpi=300, transparent=True)
    plt.show()

    # Plot b2
    plt.figure(figsize=(1, 4))
    plt.imshow(np.transpose([b2]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\b2_jgp.png",
                    dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\b2_jgp.pdf",
                    dpi=300, transparent=True)
    plt.show()

    # Plot l1
    plt.figure(figsize=(1, 4))
    plt.imshow(np.transpose([l1]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\l1_jgp.png",
                    dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\l1_jgp.pdf",
                    dpi=300, transparent=True)
    plt.show()

    # Plot l2
    plt.figure(figsize=(1, 4))
    plt.imshow(np.transpose([l2]), aspect=0.1, cmap="viridis", interpolation="none")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if True:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\l2_jgp.png",
                    dpi=300, transparent=True)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\l2_jgp.pdf",
                    dpi=300, transparent=True)
    plt.show()

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    rcParams['ps.useafm']=True
    rcParams['ps.fonttype']=42
    rcParams['pdf.fonttype']=42
    rcParams['font.family'] = 'Arial'

    xtick_size = 24
    ytick_size = 24
    xlabel_size = 32
    ylabel_size = 24

    # current_mass = 2000
    # current_charge = 5
    current_mass = 5978
    current_charge = 7
    # current_mass = 1355
    # current_charge = 1
    current_phaseres = 4

    mock_data = generate_mock_data(current_mass, current_charge)
    plot_wave_under_data(mock_data, current_charge, phaseres=current_phaseres)
    emat = get_full_emat(mock_data, phaseres=current_phaseres)
    plot_emat_highlight_z(emat, current_charge, phaseres=current_phaseres)
    # # predict_plot_predvec(mock_data, phaseres=current_phaseres)
    # # plot_nn_layers(emat)




