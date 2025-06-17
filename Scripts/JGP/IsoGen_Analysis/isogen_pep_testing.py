from unidec.IsoDec.IsoGen.isogenc import *
from unidec.IsoDec.IsoGen.isogenmass import IsoGenMassEngine
from unidec.IsoDec.IsoGen.isogenpep import *
from unidec.IsoDec.IsoGen.isogen_tools import *
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenPepEngine()
eng_avg = IsoGenMassEngine()

if False:
    masses = [1000, 5000, 10000, 50000]
    isolen = 128

    for mass in masses:
        py_fft_dist = pepmass_to_dist(mass, isolen=128)
        py_nn_dist = eng_avg.predict(mass)
        c_fft_dist = fft_pep_mass_to_dist(mass, isolen=128)
        c_nn_dist = nn_pep_mass_to_dist(mass)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.legend()
        plt.show()

if True:
    seqs = ["ADAKCGS", "AACCGKHIYHHHCCHG"]

    for seq in seqs:
        py_fft_dist = peptide_to_dist(seq, 64)
        py_nn_dist = eng.predict(seq)
        c_fft_dist = fft_pep_seq_to_dist(seq)
        c_nn_dist = nn_pep_seq_to_dist(seq)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.legend()
        plt.show()


