from unidec.IsoDec.IsoGen.isogenc import *
from matplotlib.pyplot import rcParams
from unidec.IsoDec.IsoGen.isogen_tools import rnaseq_to_mass


def make_nice_dist_plot(fft_dist, nn_seq_dist, nn_mass_dist, save_figs=False):
    tick_size = 16
    label_size = 20
    linewidth = 2

    fig, ax = plt.subplots(figsize=(6, 3))
    # for i in range(len(fft_dist)):
    #     plt.vlines(i, 0, fft_dist[i], color='k', linewidth=linewidth)
    plt.plot(fft_dist, color='k', linewidth=linewidth)
    plt.plot(nn_seq_dist, color='b', linewidth=linewidth)
    plt.plot(nn_mass_dist, color='r', linewidth=linewidth)
    plt.ylim(bottom=0)
    plt.xlim(30,130)
    plt.xticks(fontsize=tick_size)
    plt.yticks([0,0.5,1], fontsize=tick_size)
    plt.xlabel('Isotopologue Number', fontsize=label_size)
    plt.ylabel("Abundance", fontsize=label_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\polyA.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\polyA.png", dpi=300, transparent=True)
    plt.show()

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'

    tick_size = 16
    label_size = 20
    linewidth = 3

    seq1 = 'A' * 500
    seq2 = 'U' * 500

    mass = rnaseq_to_mass(seq1)

    fft_dist = fft_rna_seq_to_dist(seq1, isolen=128)
    nn_seq_dist = nn_rna_seq_to_dist(seq1, isolen=128)
    nn_mass_dist = nn_rna_mass_to_dist(mass, isolen=128)

    make_nice_dist_plot(fft_dist, nn_seq_dist, nn_mass_dist, save_figs=False)