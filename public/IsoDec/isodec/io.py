"""Spectrum readers with optional UniDec support for complex formats."""

from pathlib import Path

import numpy as np


_STANDALONE_TYPES = [".txt", ".dat", ".csv", ".npz"]
_UNIDEC_TYPES = [".raw", ".mzxml", ".mzml", ".mzml.gz", ".gz", ".i2ms", ".dmt", ".bin"]
recognized_types = _UNIDEC_TYPES + _STANDALONE_TYPES


class TextSpectrumImporter:
    """Minimal single-scan reader for NumPy-compatible text and NPZ files."""

    def __init__(self, path, **_):
        self.path = Path(path)
        self.scans = np.array([1])
        self.centroided = True
        self.polarity = "Positive"
        self._data = self._load()

    def _load(self):
        if self.path.suffix.lower() == ".npz":
            with np.load(self.path) as archive:
                if not archive.files:
                    raise ValueError(f"No arrays found in {self.path}")
                data = archive[archive.files[0]]
        else:
            delimiter = "," if self.path.suffix.lower() == ".csv" else None
            data = np.loadtxt(self.path, delimiter=delimiter)
        data = np.asarray(data, dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(f"Spectrum must contain at least two columns: {self.path}")
        return np.ascontiguousarray(data[:, :2])

    def get_single_scan(self, _scan):
        return self._data.copy()

    def check_centroided(self):
        if len(self._data) < 4:
            self.centroided = True
        else:
            spacing = np.diff(self._data[:, 0])
            median = np.median(spacing)
            self.centroided = not np.allclose(spacing, median, rtol=0.05, atol=0)
        return self.centroided

    def get_scan_time(self, _scan):
        return 0.0

    def get_ms_order(self, _scan):
        return 1


def _unidec_importer_factory():
    try:
        from unidec.UniDecImporter.ImporterFactory import ImporterFactory as UniDecImporterFactory
    except ImportError as error:
        raise ImportError(
            "This spectrum format requires UniDec. Install the optional "
            "integration with `pip install isodec[unidec]`."
        ) from error
    return UniDecImporterFactory


class ImporterFactory:
    """Create built-in text readers or delegate complex formats to UniDec."""

    @staticmethod
    def create_importer(file_path, **kwargs):
        path = Path(file_path)
        suffix = path.suffix.lower()
        if suffix in _STANDALONE_TYPES:
            return TextSpectrumImporter(path, **kwargs)
        return _unidec_importer_factory().create_importer(file_path, **kwargs)


__all__ = ["ImporterFactory", "TextSpectrumImporter", "recognized_types"]
