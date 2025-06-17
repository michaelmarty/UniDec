from unidec.IsoDec.IsoGen.isogenrna import *
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogenrna_averagine import IsoGenRNAveragineEngine

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenRNAEngine()
eng_avg = IsoGenRNAveragineEngine()


if True:
    print("Testing RNA mass to distribution methods...")
    masses = [1000, 5000, 10000, 50000, 150000]
    isolen = 128

    for mass in masses:
        py_fft_dist = rnamass_to_dist(mass, isolen=128)
        py_nn_dist = eng_avg.predict(mass)
        c_fft_dist = fft_rna_mass_to_dist(mass, isolen=128)
        c_nn_dist = nn_rna_mass_to_dist(mass)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.title(str(mass))
        plt.legend()
        plt.show()

if True:
    print("Testing RNA sequence to distribution methods...")
    testseqs = ["GUAC", "GUACGUACGUACGUAC", "GUACGUACGUACGUACGUACGUACGUACGUACGUAC", "GUACGUACGUACGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC", "AAAAAAAAAAAAA", "UUUUUUUUUUUUU"]
    isolen = 128

    for seq in testseqs:
        py_fft_dist = rnaseq_to_dist(seq, isolen=isolen)
        py_nn_dist = eng.predict(seq)
        c_fft_dist = fft_rna_seq_to_dist(seq, isolen=isolen)
        c_nn_dist = nn_rna_seq_to_dist(seq)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.legend()
        plt.show()