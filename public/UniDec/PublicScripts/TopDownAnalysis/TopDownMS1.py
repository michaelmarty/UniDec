from pyopenms import *
import unidec.engine as unidec
import time
import numpy as np
import unidec.tools as ud


def ms1_processing(path, outpath=None, fullms1=True, agressive_speed=False):
    # Start timer
    begin = time.perf_counter()
    # Set outpath
    if outpath is None:
        outpath = path[:-8] + "_deconvolved"
    # Load file
    exp = MSExperiment()
    MzMLFile().load(path, exp)

    # Check number of scans
    num_scans = exp.getNrSpectra()
    print("Number of Scans:", num_scans)

    # Set up output experiment
    out_exp = MSExperiment()

    # Create a UniDec engine
    eng = unidec.UniDec()

    # Set some defaults
    eng.config.numit = 10
    if agressive_speed:
        print("Running agressive settings")
        # Bin together every 10 datapoints to drop the resolution by 10x
        eng.config.mzbins = 10
        # Remove 80% of the data
        eng.config.reductionpercent = 80
        # Set an intensity threshold of 1e6
        eng.config.intthresh = 1e6
    # Turn off normalization
    eng.datanorm = 0
    eng.peaknorm = 0

    # Loop through spectra
    for i, spec in enumerate(exp):
        # Get MS level
        ms_level = spec.getMSLevel()
        if ms_level != 1:
            # Append to output experiment
            out_exp.addSpectrum(spec)
        else:
            print("MS Level:", ms_level, i)
            # Get the m/z and intensity values
            mz, intensity = spec.get_peaks()
            inputdata = np.transpose([mz, intensity])
            inputdata = ud.dataprep(inputdata, eng.config)

            if len(inputdata) < 2:
                continue

            # Run UniDec
            eng.pass_data_in(inputdata, silent=True, refresh=True)
            eng.process_data(silent=True)
            eng.run_unidec(silent=True, efficiency=True)

            if fullms1:
                # Get the mass and intensity values
                try:
                    mass, intensity = eng.data.massdat[:, 0], eng.data.massdat[:, 1]
                except:
                    mass, intensity = [], []

                # Update the old spectrum
                spec2 = spec
                spec2.set_peaks((mass, intensity))
                out_exp.addSpectrum(spec2)

            else:
                # Create a new spectrum
                spec2 = spec
                try:
                    peaks = eng.pick_peaks()
                    spec2.set_peaks((peaks[:, 0], peaks[:, 1]))
                except:
                    spec2.set_peaks(([], []))


                out_exp.addSpectrum(spec2)

    # Write the output experiment to a file
    MzMLFile().store(outpath + ".mzML", out_exp)


    print("Total Time:", time.perf_counter() - begin)


if __name__ == "__main__":
    # set paths
    path = "C:\\Data\\TabbData\\20170105_L_MaD_ColC4_Ecoli20161108-10_01_ms1.mzML.gz"

    # Run the function
    ms1_processing(path, fullms1=False, agressive_speed=False)
