import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from pyopenms import *
from unidec.modules import hdf5_tools as hd


def write_hdf5_to_mzml(hdf5_path, outfile=None, plot=False):
    """
    Write the HDF5 data from MetaUniDec or UniChrom to an mzML file
    :param hdf5_path: the path to the HDF5 file
    :param outfile: optional, the path to the output file. Will default to the same path as the HDF5 file with the extension changed to _deconvolved.mzML
    :param plot: whether to plot the data
    :return: None
    """
    if outfile is None:
        outfile = os.path.splitext(hdf5_path)[0] + "_deconvolved.mzML"

    # Load the HDF5 data
    massaxis = hd.get_data(hdf5_path, "ms_dataset", "mass_axis")
    massgrid = hd.get_data(hdf5_path, "ms_dataset", "mass_grid")

    # Reshape the data
    num = int(len(massgrid) / len(massaxis))
    massgrid = massgrid.reshape(int(num), len(massaxis))
    print("Data Shape:", np.shape(massgrid))

    # Get the scan start times
    try:
        times = hd.get_metadata(hdf5_path, "timestart")
    except:
        times = np.arange(0, num)


    # Write the data to an mzML file
    exp = MSExperiment()
    for i in range(0, num):
        spec = MSSpectrum()
        spec.setRT(times[i])
        spec.set_peaks((massaxis, massgrid[i]))
        exp.addSpectrum(spec)
    MzMLFile().store(outfile, exp)
    print("Wrote Data to:", outfile)

    # plot the data
    if plot:
        for d in massgrid:
            plt.plot(massaxis, d)
        plt.show()


if __name__ == "__main__":
    # Path to the HDF5 file
    hdf5_path = "../unidec/bin/Example Data/UniChrom/SEC_Native_Herceptin.hdf5"
    # Write the data to an mzML file
    write_hdf5_to_mzml(hdf5_path, plot=True)
