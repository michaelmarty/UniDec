import os
import pickle
import time

import numpy as np
from unidec.IsoDec import encoding
import unidec.tools as ud
import matplotlib.pyplot as plt
import matplotlib as mpl
import torch


def load_pkl_training_data(file_path):
    """
    Load the complete peak collection from a pickle file.

    Parameters:
        file_path (str): Path to the .pkl file.
    """
    with open(file_path, "rb") as infile_handle:
        return pickle.load(infile_handle)

def encode_phase_file_profile(file_path, onedropper=0.0, types=None):
    if types is None:
        types = [[1, 0], [4, 0], [4, 1], [8, 0], [8, 1]]

    training = []
    test = []
    zdist = []

    try:
        with open(file_path, "rb") as f:
            pks = pickle.load(f)
    except Exception as e:
        print("Error Loading File:", file_path, e)
        return training, test, zdist, types

    # Encode each centroid
    print("File:", file_path, len(pks))
    for p in pks:
        z = p.z
        if z == 1:
            # toss it out with a onedropper% chance
            r = np.random.rand()
            if r < onedropper:
                continue

        zdist.append(z)

        centroid = p.centroids
        profile = p.profile_data

        if len(centroid) < 3 or len(profile) < 3:
            print("Skipping peak with insufficient data:", len(centroid), len(profile))
            continue

        cemats = []
        pemats = []
        for t in types:
            phaseres = t[0]
            shift = t[1]

            emat = encoding.encode_phase(centroid, phaseres=phaseres, shift=shift)
            emat = torch.as_tensor(emat, dtype=torch.float32)
            cemats.append(emat)

            pemat = encoding.encode_phase(profile, phaseres=phaseres, shift=shift)
            pemat = torch.as_tensor(pemat, dtype=torch.float32)
            pemats.append(pemat)

        # Flatten emats and pemats into a single tensor for each, add profile and centroid too with z
        output_tuple = (*cemats, *pemats, centroid, profile, z)
        # randomly sort into training and test data
        r = np.random.rand()
        if r < 0.9:
            training.append(output_tuple)
        else:
            test.append(output_tuple)
    return training, test, zdist, types

def save_encoding(outdata, outfile, types):
    num_types = len(types)
    expected_record_length = 2 * num_types + 3
    if any(len(record) != expected_record_length for record in outdata):
        raise ValueError(
            f"Each encoded record must contain {expected_record_length} values for {num_types} encoding types"
        )

    output_data = {}
    for i, (phaseres, shift) in enumerate(types):
        shift = int(bool(shift))
        centroid_key = f"emat_c_{phaseres}_{shift}"
        profile_key = f"emat_p_{phaseres}_{shift}"
        if centroid_key in output_data or profile_key in output_data:
            raise ValueError(f"Duplicate encoding type: phaseres={phaseres}, shift={shift}")

        if outdata:
            output_data[centroid_key] = torch.stack([record[i] for record in outdata]).cpu().numpy()
            output_data[profile_key] = torch.stack(
                [record[num_types + i] for record in outdata]
            ).cpu().numpy()
        else:
            output_data[centroid_key] = np.empty((0, 50, phaseres), dtype=np.float32)
            output_data[profile_key] = np.empty((0, 50, phaseres), dtype=np.float32)

    centroids = np.empty(len(outdata), dtype=object)
    profiles = np.empty(len(outdata), dtype=object)
    centroids[:] = [record[-3] for record in outdata]
    profiles[:] = [record[-2] for record in outdata]
    output_data["centroids"] = centroids
    output_data["profiles"] = profiles
    output_data["z"] = np.fromiter((record[-1] for record in outdata), dtype=np.int32, count=len(outdata))

    print("Saving to:", outfile, "Length:", len(outdata))
    np.savez_compressed(outfile, **output_data)


def encode_dir(pkldir, outdir=None, name="profile", file_string="_profile_data.pkl", maxfiles=None, plot=False, types=None, **kwargs):
    starttime = time.perf_counter()
    training = []
    test = []
    zdist = []

    if types is None:
        types = [[1, 0], [4, 0], [4, 1], [8, 0], [8, 1]]
    encoding_types = [list(encoding_type) for encoding_type in types]

    files = ud.match_files_recursive(pkldir, file_string)
    if maxfiles is not None:
        files = files[:maxfiles]

    # print("Files:", files)
    for file in files:
        tr, te, zd, file_types = encode_phase_file_profile(file, types=encoding_types, **kwargs)
        if [list(encoding_type) for encoding_type in file_types] != encoding_types:
            raise ValueError(f"Encoding types changed while processing {file}")

        training.extend(tr)
        test.extend(te)
        zdist.extend(zd)


    output_directory = outdir if outdir is not None else pkldir
    os.makedirs(output_directory, exist_ok=True)
    training_file = os.path.join(output_directory, f"training_data_{name}.npz")
    test_file = os.path.join(output_directory, f"test_data_{name}.npz")
    save_encoding(training, training_file, encoding_types)
    save_encoding(test, test_file, encoding_types)

    elapsed = time.perf_counter() - starttime
    print("Time:", elapsed, "Number of peaks:", len(zdist))
    if plot:
        plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
        plt.show()

    return training_file, test_file, zdist


def encode_dir_list(topdir, dirlist, name="profile", plot=True, check_existing=True, **kwargs):
    os.chdir(topdir)
    zdist = []
    for subdir in dirlist:
        full_path = os.path.join(topdir, subdir)
        if check_existing:
            training_file = os.path.join(full_path, f"training_data_{name}.npz")
            test_file = os.path.join(full_path, f"test_data_{name}.npz")
            if os.path.exists(training_file) and os.path.exists(test_file):
                print(f"Skipping {subdir} (already processed)")
                continue
        tr_file, te_file, zd = encode_dir(full_path, name=name, plot=False, **kwargs)
        zdist.extend(zd)

    print("Total Number of Peaks:", len(zdist))
    if plot:
        plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
        plt.savefig(os.path.join(topdir, f"zdist_{name}.png"))
        plt.savefig(os.path.join(topdir, f"zdist_{name}.pdf"))
        plt.show()

    return zdist

data_dirs = ["MSV000090488",
             "MSV000091923",
             "PXD045560",
             "PXD046651",
             "PXD027650",
             "PXD041357",
             "PXD042921",
             "PXD069439"
             ]



if __name__ == "__main__":

    mpl.use("WxAgg")

    topdir = r"Z:\Group Share\JGP"
    encode_dir_list(topdir, data_dirs, check_existing=True)
    exit()

    topdir = r"Z:\Group Share\JGP\PXD069439"
    encode_dir(topdir, name="profile", plot=True)
    exit()


    file_path=r"Z:\Group Share\JGP\PXD069439\TDP_HH_F12_profile_data.pkl"

    training, test, zdist, types = encode_phase_file_profile(file_path, onedropper=0.1)

    plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
    plt.show()

    exit()

    # Load the profile data
    pks = load_pkl_training_data(file_path)

    # Plot the first peak in both centroid and profile
    for p in pks:
        if p.z > 3:
            first_peak = p
            # Replace centroids with a pure gaussian peak centered at 0 with a width of 0.01
            #
            # sigma = 0.1
            # x = np.linspace(-3*sigma, 3*sigma, 1000)
            # y = np.exp(-((x) ** 2) / (2 * (sigma ** 2)))
            # y /= np.max(y)  # Normalize to max of 1
            # first_peak.centroids = np.column_stack((x, y))

            plt.figure(figsize=(10, 5))
            plt.subplot(1, 3, 1)
            plt.plot(first_peak.centroids[:, 0], first_peak.centroids[:, 1], "o", label='Centroid Data')
            plt.plot(first_peak.profile_data[:, 0], first_peak.profile_data[:, 1], label='Profile Data')

            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.legend()


            phaseres=8
            shift=True
            # For first peak, plot the encoded centroid data in the second panel using an imshow plot
            plt.subplot(1, 3, 2)
            emat = encoding.encode_phase(first_peak.centroids, phaseres=phaseres, shift=shift)
            plt.imshow(emat, cmap='viridis', aspect='auto')

            # For the last panel, repeat but with profile data
            plt.subplot(1, 3, 3)
            emat_profile = encoding.encode_phase(first_peak.profile_data, phaseres=phaseres, shift=shift)
            plt.imshow(emat_profile, cmap='viridis', aspect='auto')

            plt.title(str(first_peak.z))
            plt.show()
