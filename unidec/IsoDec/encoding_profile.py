import pickle

from unidec.IsoDec import encoding
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("WxAgg")

def load_pkl_training_data(file_path):
    """
    Load the complete peak collection from a pickle file.

    Parameters:
        file_path (str): Path to the .pkl file.
    """
    with open(file_path, "rb") as infile_handle:
        return pickle.load(infile_handle)

if __name__ == "__main__":
    file_path=r"Z:\Group Share\JGP\PXD069439\TDP_HH_F12_profile_data.pkl"

    # Load the profile data
    pks = load_pkl_training_data(file_path)

    # Plot the first peak in both centroid and profile
    for p in pks:
        if p.z > 3:
            first_peak = p

            plt.figure(figsize=(10, 5))
            plt.subplot(1, 3, 1)
            plt.plot(first_peak.centroids[:, 0], first_peak.centroids[:, 1], "o", label='Centroid Data')
            plt.plot(first_peak.profile_data[:, 0], first_peak.profile_data[:, 1], label='Profile Data')

            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.legend()


            phaseres=8
            # For first peak, plot the encoded centroid data in the second panel using an imshow plot
            plt.subplot(1, 3, 2)
            emat = encoding.encode_phase(first_peak.centroids, phaseres=phaseres, shift=True)
            plt.imshow(emat, cmap='viridis', aspect='auto')

            # For the last panel, repeat but with profile data
            plt.subplot(1, 3, 3)
            emat_profile = encoding.encode_phase(first_peak.profile_data, phaseres=phaseres, shift=True)
            plt.imshow(emat_profile, cmap='viridis', aspect='auto')

            plt.title(str(first_peak.z))
            plt.show()
