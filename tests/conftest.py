from pathlib import Path

import numpy as np
import pytest


TEST_SPECTRUM = Path(__file__).with_name("test_spectrum.txt")


@pytest.fixture(scope="session")
def spectrum():
    data = np.loadtxt(TEST_SPECTRUM)
    data.setflags(write=False)
    return data


@pytest.fixture(scope="session")
def processed_spectrum(spectrum):
    from isodec import IsoDecRuntime

    runtime = IsoDecRuntime()
    return runtime.batch_process_spectrum(spectrum, centroided=True, refresh=True)
