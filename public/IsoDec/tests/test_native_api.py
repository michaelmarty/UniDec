import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from isodec import IsoDecConfig, IsoDecWrapper, __version__
from isodec.c_interface import default_dll_path, example


@pytest.mark.parametrize(
    "module_name",
    [
        "__main__.py",
        "altdecon.py",
        "c_interface.py",
        "match.py",
        "plots.py",
        "runtime.py",
        "trainingdata_profile.py",
    ],
)
def test_core_runnable_modules_load_without_package_context(tmp_path, module_name):
    module_path = Path(__file__).parents[1] / "isodec" / module_name
    probe = (
        "import runpy, sys; "
        "runpy.run_path(sys.argv[1], run_name='direct_execution_test')"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe, str(module_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_all_runnable_modules_with_relative_imports_bootstrap_the_package():
    package_dir = Path(__file__).parents[1] / "isodec"
    for module_path in package_dir.glob("*.py"):
        source = module_path.read_text(encoding="utf-8")
        is_runnable = 'if __name__ == "__main__"' in source
        has_relative_imports = any(line.startswith("from .") for line in source.splitlines())
        if is_runnable and has_relative_imports:
            assert "if not __package__:" in source, module_path.name
            assert '__package__ = "isodec"' in source, module_path.name


def test_native_library_is_packaged():
    assert Path(default_dll_path).is_file()


def test_default_configuration_is_standalone_and_valid():
    config = IsoDecConfig()
    assert config.phaseres == 8
    assert config.minpeaks == 3
    assert config.adductmass == pytest.approx(1.007276467)
    assert config.mzwindowlb < 0 < config.mzwindowub


def test_native_example_charge_and_deconvolution():
    wrapper = IsoDecWrapper()
    assert wrapper.predict_charge(example) == 11

    peaks = wrapper.process_spectrum(example)
    assert [peak.z for peak in peaks] == [11, 5]
    np.testing.assert_allclose(
        [peak.monoiso for peak in peaks],
        [6250.5923, 2832.5012],
        rtol=0,
        atol=0.02,
    )
    assert all(peak.matchedintensity > 0 for peak in peaks)
