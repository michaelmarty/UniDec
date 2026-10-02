import ctypes
import time
from importlib import metadata
import os
import platform
import sys
from pathlib import Path

import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "isodec"

from .match import MatchedPeak, MatchedCollection, MatchedMass
from .config import IsoDecConfig
from .plots import cplot

_system = platform.system()


def _group_array(value, name, pairs=True, allow_empty=False):
    """Validate lengths and shapes before exposing a NumPy buffer to C."""
    array = np.ascontiguousarray(value, dtype=np.float64)
    if (array.ndim != (2 if pairs else 1) or
            (pairs and array.shape[1] != 2)):
        raise ValueError(f"{name} must have shape (N, 2)" if pairs
                         else f"{name} must be one-dimensional")
    if (not allow_empty and len(array) == 0) or array.size > np.iinfo(np.int32).max:
        raise ValueError(f"{name} has an invalid size")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array

_library_names = {
    "Windows": "isodeclib.dll",
    "Linux": "isodeclib.so",
    "Darwin": "isodeclib.dylib",
}
try:
    dllname = _library_names[_system]
except KeyError as error:
    raise ImportError(f"IsoDec does not support {_system!r}") from error

_package_dir = Path(__file__).resolve().parent
_packaged_library_path = _package_dir / "bin" / dllname
_source_library_path = _package_dir / dllname


def _editable_library_path():
    """Locate scikit-build's native install during an editable install."""
    try:
        candidate = Path(
            metadata.distribution("isodec").locate_file(
                Path("isodec") / "bin" / dllname
            )
        )
    except metadata.PackageNotFoundError:
        return None
    return candidate if candidate.is_file() else None

if _packaged_library_path.is_file():
    default_dll_path = _packaged_library_path
elif (_package_dir.parent / "pyproject.toml").is_file() and _source_library_path.is_file():
    default_dll_path = _source_library_path
elif _editable_library_path() is not None:
    default_dll_path = _editable_library_path()
else:
    default_dll_path = _packaged_library_path

example = np.array(
    [
        [5.66785531e02, 1.47770838e06],
        [5.67057354e02, 1.54980838e06],
        [5.67507468e02, 5.21600520e07],
        [5.67708173e02, 8.35557760e07],
        [5.67908401e02, 7.28264240e07],
        [5.68060254e02, 1.87337225e06],
        [5.68108674e02, 4.35435520e07],
        [5.68239256e02, 3.88155375e06],
        [5.68309390e02, 2.05468060e07],
        [5.68509951e02, 7.18109250e06],
        [5.68707871e02, 2.30373500e06],
        [5.69150563e02, 1.57598062e06],
        [5.69243121e02, 1.96390440e07],
        [5.69334393e02, 6.82677120e07],
        [5.69425337e02, 1.22867432e08],
        [5.69516492e02, 1.45702336e08],
        [5.69607541e02, 1.20801936e08],
        [5.69698595e02, 1.06786072e08],
        [5.69789906e02, 6.56232960e07],
        [5.69881208e02, 3.41013880e07],
        [5.69972168e02, 1.70930360e07],
        [5.70063432e02, 9.17621100e06],
        [5.70699369e02, 1.96462650e06],
    ]
)

isodist = ctypes.c_float * 64
matchedinds = ctypes.c_int * 32


# print(isodist)
class MPStruct(ctypes.Structure):
    _fields_ = [
        ("mz", ctypes.c_float),
        ("z", ctypes.c_int),
        ("monoiso", ctypes.c_float),
        ("peakmass", ctypes.c_float),
        ("avgmass", ctypes.c_float),
        ("area", ctypes.c_float),
        ("peakint", ctypes.c_float),
        ("matchedindsiso", ctypes.c_int * 64),
        ("matchedindsexp", ctypes.c_int * 64),
        ("isomz", isodist),
        ("isodist", isodist),
        ("isomass", isodist),
        ("monoisos", ctypes.c_float * 16),
        ("startindex", ctypes.c_int),
        ("endindex", ctypes.c_int),
        ("score", ctypes.c_float),
        ("realisolength", ctypes.c_int),
    ]


class FragmentHitStruct(ctypes.Structure):
    _fields_ = [
        ("fragment_index", ctypes.c_int),
        ("charge", ctypes.c_int),
        ("left", ctypes.c_int),
        ("right", ctypes.c_int),
        ("isotope_count", ctypes.c_int),
        ("match_count", ctypes.c_int),
        ("score", ctypes.c_double),
        ("scale", ctypes.c_double),
        ("isotope_indexes", ctypes.c_int * 128),
        ("centroid_indexes", ctypes.c_int * 128),
    ]


class MassGroupPeakStruct(ctypes.Structure):
    _fields_ = [
        ("monoiso", ctypes.c_double), ("mz", ctypes.c_double),
        ("peakint", ctypes.c_double), ("matchedintensity", ctypes.c_double),
        ("avgmass", ctypes.c_double), ("rt", ctypes.c_double),
        ("charge", ctypes.c_int), ("scan", ctypes.c_int),
        ("monoisos", ctypes.POINTER(ctypes.c_double)), ("monoisos_count", ctypes.c_int),
        ("massdist", ctypes.POINTER(ctypes.c_double)), ("massdist_count", ctypes.c_int),
        ("decon_centroids", ctypes.POINTER(ctypes.c_double)), ("centroid_count", ctypes.c_int),
        ("float32_intensity", ctypes.c_int), ("float32_massdist", ctypes.c_int),
    ]


class MassGroupResultStruct(ctypes.Structure):
    _fields_ = [
        ("monoiso", ctypes.c_double), ("lookup_mass", ctypes.c_double),
        ("apexintensity", ctypes.c_double),
        ("totalintensity", ctypes.c_double), ("avgmass", ctypes.c_double),
        ("minrt", ctypes.c_double), ("maxrt", ctypes.c_double),
        ("apexrt", ctypes.c_double), ("minscan", ctypes.c_int),
        ("maxscan", ctypes.c_int), ("apexscan", ctypes.c_int),
        ("seed_index", ctypes.c_int), ("totalpeaks", ctypes.c_int),
        ("monoisos", ctypes.POINTER(ctypes.c_double)), ("monoisos_count", ctypes.c_int),
        ("massdist", ctypes.POINTER(ctypes.c_double)), ("massdist_count", ctypes.c_int),
        ("decon_centroids", ctypes.POINTER(ctypes.c_double)), ("centroid_count", ctypes.c_int),
    ]


class IDSettings(ctypes.Structure):
    _fields_ = [
        ("phaseres", ctypes.c_int),
        ("verbose", ctypes.c_int),
        ("peakwindow", ctypes.c_int),
        ("peakthresh", ctypes.c_float),
        ("minpeaks", ctypes.c_int),
        ("css_thresh", ctypes.c_float),
        ("matchtol", ctypes.c_float),
        ("maxshift", ctypes.c_int),
        ("mzwindow", ctypes.c_float * 2),
        ("plusoneintwindow", ctypes.c_float * 2),
        ("knockdown_rounds", ctypes.c_int),
        ("min_score_diff", ctypes.c_float),
        ("minareacovered", ctypes.c_float),
        ("isolength", ctypes.c_int),
        ("mass_diff_c", ctypes.c_double),
        ("adductmass", ctypes.c_float),
        ("minusoneaszero", ctypes.c_int),
        ("isotopethreshold", ctypes.c_float),
        ("datathreshold", ctypes.c_float),
        ("zscore_threshold", ctypes.c_float),
    ]


class IDConfig(ctypes.Structure):
    _fields_ = [
        ("verbose", ctypes.c_int),
        ("pres", ctypes.c_int),
        ("maxz", ctypes.c_int),
        ("elen", ctypes.c_int),
        ("l1", ctypes.c_int),
        ("l2", ctypes.c_int),
        ("l3", ctypes.c_int),
        ("l4", ctypes.c_int),
        ("dlen", ctypes.c_int),
    ]


def config_to_settings(config):
    # print(config)
    settings = IDSettings()
    settings.phaseres = int(config.phaseres)
    settings.verbose = int(config.verbose)
    settings.peakwindow = int(config.peakwindow)
    settings.peakthresh = float(config.peakthresh)
    settings.minpeaks = int(config.minpeaks)
    settings.css_thresh = float(config.css_thresh)
    settings.matchtol = float(config.matchtol)
    settings.maxshift = int(config.maxshift)
    settings.mzwindow = (config.mzwindowlb, config.mzwindowub)
    settings.plusoneintwindow = (config.plusoneintwindowlb, config.plusoneintwindowub)
    settings.knockdown_rounds = int(config.knockdown_rounds)
    settings.min_score_diff = config.min_score_diff
    settings.minareacovered = config.minareacovered
    settings.isolength = 64
    settings.mass_diff_c = config.mass_diff_c
    settings.adductmass = config.adductmass
    settings.minusoneaszero = config.minusoneaszero
    settings.isotopethreshold = config.isotopethreshold
    settings.datathreshold = config.datathreshold
    settings.zscore_threshold = config.zscore_threshold
    return settings


class IsoDecWrapper:
    def __init__(self, dllpath=None):
        if dllpath is None:
            dllpath = default_dll_path

        dllpath = Path(dllpath).resolve()
        if not dllpath.is_file():
            raise ImportError(
                f"IsoDec's native library is missing: {dllpath}. Reinstall "
                "isodec using a compatible wheel, or build from source "
                "with CMake and a native compiler."
            )

        modelpath = _package_dir / "modelparams"

        self.modeldir = modelpath

        self._dll_directory_handle = None
        if _system == "Windows":
            self._dll_directory_handle = os.add_dll_directory(str(dllpath.parent))
        try:
            self.c_lib = ctypes.CDLL(str(dllpath))
        except OSError as error:
            raise ImportError(f"Unable to load IsoDec's native library {dllpath}: {error}") from error

        self.c_lib.encode.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float),
            IDConfig,
            IDSettings,
        ]
        self.c_lib.encode.restype = ctypes.c_int

        self.c_lib.predict_charge.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_char_p,  # Model path
        ]
        self.c_lib.predict_charge.restype = ctypes.c_int

        self.c_lib.process_spectrum.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(MPStruct),
            IDSettings,
            ctypes.c_char_p,
        ]
        self.c_lib.process_spectrum.restype = ctypes.c_int

        self.c_lib.DefaultSettings.argtypes = []
        self.c_lib.DefaultSettings.restype = IDSettings

        self._fragment_match = getattr(self.c_lib, "match_fragment_batch", None)
        if self._fragment_match is not None:
            self._fragment_match.argtypes = [
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_float),
                ctypes.c_int, ctypes.c_int, ctypes.c_int,
                ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int,
                ctypes.c_double, ctypes.c_double, ctypes.c_int,
                ctypes.POINTER(ctypes.POINTER(FragmentHitStruct)), ctypes.POINTER(ctypes.c_int),
            ]
            self._fragment_match.restype = ctypes.c_int
            self.c_lib.free_fragment_hits.argtypes = [ctypes.POINTER(FragmentHitStruct)]
            self.c_lib.free_fragment_hits.restype = None

        self._group_centroid_match = getattr(self.c_lib, "match_group_centroids_batch", None)
        if self._group_centroid_match is not None:
            self._group_centroid_match.argtypes = [
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_double,
                ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.c_int,
            ]
            self._group_centroid_match.restype = ctypes.c_int

        self._mass_group = getattr(self.c_lib, "group_mass_peaks_batch", None)
        if self._mass_group is not None:
            version = getattr(self.c_lib, "mass_group_abi_version", None)
            if version is None:
                raise ImportError("native mass-grouping ABI has no version")
            version.argtypes = []
            version.restype = ctypes.c_int
            if version() != 1:
                raise ImportError("unsupported native mass-grouping ABI version")
            self._mass_group.argtypes = [
                ctypes.POINTER(MassGroupPeakStruct), ctypes.c_int,
                ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                ctypes.c_int, ctypes.c_int,
                ctypes.POINTER(ctypes.POINTER(MassGroupResultStruct)),
                ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
            ]
            self._mass_group.restype = ctypes.c_int
            self.c_lib.free_mass_groups.argtypes = [
                ctypes.POINTER(MassGroupResultStruct), ctypes.c_int,
                ctypes.POINTER(ctypes.c_int),
            ]
            self.c_lib.free_mass_groups.restype = None

        self.modeldir = str(modelpath)
        # self.modelpath = ctypes.c_char_p(
        #     os.path.join(self.modeldir, "phase_model_8.bin").encode()
        # )
        # Create null pointer for model path
        self.modelpath = None
        self.config = IsoDecConfig()
        # self.determine_model()

    def match_fragment_batch(self, spectrum, masses, intensities, config, max_charge=None):
        """Return native fragment hit metadata, or None for an older library."""
        if self._fragment_match is None:
            return None
        mz = np.ascontiguousarray(spectrum[:, 0], dtype=np.float64)
        observed = np.ascontiguousarray(spectrum[:, 1], dtype=np.float64)
        masses = np.ascontiguousarray(masses, dtype=np.float64)
        intensities = np.ascontiguousarray(intensities, dtype=np.float32)
        if intensities.ndim != 2 or intensities.shape[0] != len(masses):
            raise ValueError("fragment masses and intensities must be aligned")
        pointer = ctypes.POINTER(FragmentHitStruct)()
        count = ctypes.c_int()
        result = self._fragment_match(
            mz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            observed.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(mz),
            masses.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            intensities.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(masses), intensities.shape[1], max_charge or 0,
            config.adductmass, config.isotopethreshold, config.matchtol,
            config.minpeaks, config.css_thresh, config.minareacovered,
            int(config.minusoneaszero), ctypes.byref(pointer), ctypes.byref(count),
        )
        if result != 0:
            raise ValueError("native fragment matching failed")
        try:
            return [(
                pointer[i].fragment_index, pointer[i].charge, pointer[i].left,
                pointer[i].right, pointer[i].score, pointer[i].scale,
                list(pointer[i].centroid_indexes[:pointer[i].match_count]),
                list(pointer[i].isotope_indexes[:pointer[i].match_count]),
            ) for i in range(count.value)]
        finally:
            self.c_lib.free_fragment_hits(pointer)

    def match_group_centroids_batch(self, spectrum, peaks, tolerance):
        """Return accepted local centroid indexes, or None for an older library."""
        if self._group_centroid_match is None or not peaks:
            return None
        spectrum = _group_array(spectrum, "spectrum")
        if len(peaks) > np.iinfo(np.int32).max // 2:
            raise ValueError("too many peaks")
        envelopes = [_group_array(p.isodist, "isodist", allow_empty=True) for p in peaks]
        total = sum(len(a) for a in envelopes)
        if total > np.iinfo(np.int32).max:
            raise ValueError("too many isotopes")
        for peak in peaks:
            if not 0 <= peak.startindex <= peak.endindex <= len(spectrum):
                raise ValueError("invalid centroid window")
        mz = np.ascontiguousarray(spectrum[:, 0], dtype=np.float64)
        intensity = np.ascontiguousarray(spectrum[:, 1], dtype=np.float64)
        windows = np.ascontiguousarray(
            [(peak.startindex, min(peak.endindex + 1, len(spectrum))) for peak in peaks],
            dtype=np.int32)
        offsets = np.zeros(len(peaks) + 1, dtype=np.int32)
        offsets[1:] = np.cumsum([len(a) for a in envelopes])
        isotope_mz = np.ascontiguousarray(
            np.concatenate([a[:, 0] for a in envelopes]), dtype=np.float64)
        counts = np.zeros(len(peaks), dtype=np.int32)
        indexes = np.empty(len(isotope_mz), dtype=np.int32)
        status = self._group_centroid_match(
            mz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            intensity.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(mz),
            windows.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            isotope_mz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(peaks),
            tolerance, counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            indexes.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(indexes),
        )
        if status != 0:
            raise ValueError("native mass-group centroid matching failed")
        ends = np.cumsum(counts)
        starts = np.concatenate(([0], ends[:-1]))
        return [indexes[start:end] for start, end in zip(starts, ends)]

    def group_mass_peaks_batch(self, peaks, config, order="original"):
        """Build MatchedMass objects in one native call, if supported."""
        if self._mass_group is None:
            return None
        if order not in ("original", "matched_intensity"):
            raise ValueError("mass grouping order must be original or matched_intensity")
        if not peaks:
            return [], np.array([])
        if len(peaks) > np.iinfo(np.int32).max // 2:
            raise ValueError("too many peaks")
        keepalive = []
        if not isinstance(config.maxshift, (int, np.integer)) or not 0 <= config.maxshift <= 64:
            raise ValueError("maxshift must be an integer between 0 and 64")
        inputs = (MassGroupPeakStruct * len(peaks))()
        for i, peak in enumerate(peaks):
            monoisos = _group_array(peak.monoisos, "monoisos", pairs=False)
            massdist = _group_array(peak.massdist, "massdist")
            for name in ("centroids", "matchedcentroids", "isodist"):
                value = getattr(peak, name)
                if value is not None:
                    _group_array(value, name, allow_empty=True)
            for name in ("z", "scan"):
                value = getattr(peak, name)
                if (not np.isfinite(value) or int(value) != value or
                        not np.iinfo(np.int32).min <= value <= np.iinfo(np.int32).max):
                    raise ValueError(f"{name} is outside the native integer range")
            # Use the same recovery and cache-refresh rules as MatchedMass.
            decon = _group_array(peak.calc_mass_dists(), "decon_centroids")
            keepalive.append((monoisos, massdist, decon))
            inputs[i] = MassGroupPeakStruct(
                peak.monoiso, peak.mz, peak.peakint, peak.matchedintensity,
                peak.avgmass, peak.rt, int(peak.z), int(peak.scan),
                monoisos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(monoisos),
                massdist.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(massdist),
                decon.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(decon),
                0, 0,
            )
        results = ctypes.POINTER(MassGroupResultStruct)()
        count = ctypes.c_int()
        ids = ctypes.POINTER(ctypes.c_int)()
        status = self._mass_group(
            inputs, len(peaks), config.matchtol, config.maxshift, config.mass_diff_c,
            config.css_thresh, 100, int(order == "matched_intensity"),
            ctypes.byref(results), ctypes.byref(count), ctypes.byref(ids),
        )
        if status != 0:
            raise ValueError("native mass grouping failed")
        try:
            assignment = np.ctypeslib.as_array(ids, shape=(len(peaks),)).copy()
            insertion = list(range(len(peaks)))
            if order == "matched_intensity":
                insertion.sort(key=lambda i: -peaks[i].matchedintensity)
            members_by_group = [[] for _ in range(count.value)]
            for index in insertion:
                members_by_group[assignment[index]].append(peaks[index])
            masses = []
            lookup = np.empty(count.value)
            for group_id in range(count.value):
                result = results[group_id]
                members = members_by_group[group_id]
                group = MatchedMass.__new__(MatchedMass)
                group.monoiso = result.monoiso
                lookup[group_id] = result.lookup_mass
                group.monoisos = np.ctypeslib.as_array(
                    result.monoisos, shape=(result.monoisos_count,)).copy()
                group.massdist = np.ctypeslib.as_array(
                    result.massdist, shape=(result.massdist_count * 2,)).copy().reshape(-1, 2)
                group.decon_centroids = np.ctypeslib.as_array(
                    result.decon_centroids, shape=(result.centroid_count * 2,)).copy().reshape(-1, 2)
                group.clusters = members
                group.totalpeaks = result.totalpeaks
                group.totalintensity = result.totalintensity
                group.apexintensity = result.apexintensity
                group.avgmass = result.avgmass
                group.minscan, group.maxscan = result.minscan, result.maxscan
                group.minrt, group.maxrt = result.minrt, result.maxrt
                group.apexscan, group.apexrt = result.apexscan, result.apexrt
                group.scans = np.array(list(dict.fromkeys(p.scan for p in members)))
                group.scan_intensities = {}
                for member in members:
                    if member.scan in group.scan_intensities:
                        group.scan_intensities[member.scan] += float(member.matchedintensity)
                    else:
                        group.scan_intensities[member.scan] = float(member.matchedintensity)
                first_by_charge = []
                seen_charges = set()
                for member in members:
                    if member.z not in seen_charges:
                        first_by_charge.append(member)
                        seen_charges.add(member.z)
                group.zs = np.array([p.z for p in first_by_charge])
                group.mzs = np.array([p.mz for p in first_by_charge])
                group.mzints = np.array([p.peakint for p in first_by_charge])
                distributions = [p.isodist for p in first_by_charge if p.isodist is not None]
                group.isodists = np.vstack(distributions) if distributions else None
                masses.append(group)
            return masses, lookup
        finally:
            self.c_lib.free_mass_groups(results, count, ids)

    def encode(self, centroids, maxz=50, phaseres=8, config=None):
        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        elen = maxz * phaseres
        emat = np.zeros(elen).astype(np.float32)
        idconf = IDConfig()
        idconf.pres = phaseres
        idconf.maxz = maxz
        idconf.elen = elen

        if config is not None:
            self.config = config
        settings = config_to_settings(self.config)

        self.c_lib.encode(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(len(cmz)),
            emat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            idconf,
            settings,
        )
        # Convert emat to numpy
        emat = np.ctypeslib.as_array(emat)
        return emat

    def predict_charge(self, centroids, config=None):

        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        # charge = ctypes.c_int(0)
        # self.modelpath = ctypes.c_char_p(None)

        charge = self.c_lib.predict_charge(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(len(cmz)),
            self.modelpath,
            # ctypes.byref(charge),
        )
        return charge

    def process_spectrum(self, centroids, pks=None, config=None, input_type=None):

        if input_type is None:
            input_type = "Peptide"
        type_c = ctypes.c_char_p(input_type.encode('utf-8'))

        cmz = centroids[:, 0].astype(np.double)
        cint = centroids[:, 1].astype(np.float32)
        n = len(cmz)

        if config is not None:
            self.config = config
        settings = config_to_settings(self.config)

        elems = (MPStruct * n)()
        matchedpeaks = ctypes.cast(elems, ctypes.POINTER(MPStruct))

        t1 = time.perf_counter()
        nmatched = self.c_lib.process_spectrum(
            cmz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cint.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(n),
            self.modelpath,
            matchedpeaks,
            settings,
            type_c,
        )
        t2 = time.perf_counter()
        print(f"Processed {n} centroids in {t2 - t1:.4f} seconds, found {nmatched} matched peaks.")

        if pks is None:
            pks = MatchedCollection()
        pending = []
        for p in matchedpeaks[:nmatched]:
            if p.z == 0:
                continue

            # Extract peak centroid data
            startindex = p.startindex
            endindex = p.endindex
            if startindex < 0 or endindex < 0 or endindex <= startindex:
                raise ValueError("Invalid start or end index in matched peak.")
            pkcent = centroids[startindex:endindex + 1]

            pk = MatchedPeak(p.z, p.mz, p.avgmass, centroids=pkcent, config=self.config)
            pk.monoiso = p.monoiso
            pk.peakmass = p.peakmass
            pk.avgmass = p.avgmass

            monoisos = np.array(p.monoisos)
            monoisos = monoisos[monoisos > 0]
            pk.monoisos = monoisos

            if config is not None:
                pk.scan = config.activescan
                pk.ms_order = config.activescanorder
                pk.rt = config.activescanrt

            if not 0 < p.realisolength <= len(p.isodist):
                raise ValueError("invalid native isotope length")
            isodist = np.array(p.isodist[:p.realisolength], dtype=np.float32)
            isomz = np.array(p.isomz[:p.realisolength], dtype=np.float32)
            isomass = np.array(p.isomass[:p.realisolength], dtype=np.float32)
            b1 = isodist > np.amax(isodist) * 0.001
            isodist = isodist[b1]
            isomz = isomz[b1]
            isomass = isomass[b1]
            pk.matchedintensity = np.sum(isodist)
            pk.peakint = p.peakint

            # pk.isodist = calc_isotope_dist(p.monoiso, p.z)
            pk.isodist = np.transpose((isomz, isodist))
            # pk.isodist[:, 1] = p.peakint

            pk.massdist = np.transpose((isomass, isodist))

            pk.startindex = p.startindex
            pk.endindex = p.endindex

            pending.append(pk)
        matched = self.match_group_centroids_batch(centroids, pending, self.config.matchtol)
        if matched is not None:
            for peak, indexes in zip(pending, matched):
                if len(indexes):
                    peak.matchedcentroids = peak.centroids[np.unique(indexes)]
        order = getattr(self.config, "mass_group_order", "original")
        if order not in ("original", "matched_intensity"):
            raise ValueError("mass grouping order must be original or matched_intensity")
        grouped = (self.group_mass_peaks_batch(pending, self.config, order)
                   if not pks.peaks else None)
        for pk in pending:
            pks.add_peak(pk)
        if grouped is not None:
            pks.masses, pks.monoisos = grouped
        else:
            insertion = (sorted(pending, key=lambda p: -p.matchedintensity)
                         if order == "matched_intensity" else pending)
            for pk in insertion:
                pks.add_pk_to_masses(pk, config=self.config)
        return pks

    def determine_model(self, default=True):
        if default:
            self.modelpath = None
        else:
            print("Model Directory:", self.modeldir)
            if self.config.phaseres == 4:
                self.modelpath = ctypes.c_char_p(
                    os.path.join(self.modeldir, "phase_model_4.bin").encode('utf-8')
                )
            elif self.config.phaseres == 8:

                self.modelpath = ctypes.c_char_p(
                    os.path.join(self.modeldir, "phase_model_8.bin").encode('utf-8')
                )
            else:
                print("Invalid phase resolution.", self.config.phaseres)
                raise ValueError("Invalid phase resolution.")

            if self.config.verbose:
                print(
                    "Running C code with phaseres of:",
                    self.config.phaseres,
                )


if __name__ == "__main__":
    eng = IsoDecWrapper()
    #
    # eng.config.phaseres = 4
    # eng.encode(example)
    # print(eng.predict_charge(example))
    # dat = eng.process_spectrum(example)
    #
    # print(dat)
    filepath = "C:\\Data\\IsoNN\\test2.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    pks = eng.process_spectrum(spectrum)
    print(len(pks.peaks))

    # Find peak with monoiso near 6396
    for pk in pks.masses:
        if 6395 < pk.monoiso < 6397:
            print("Found Peak:", pk)
            cplot(pk.decon_centroids, color="k", factor=1)
            cplot(pk.massdist, color="k", factor=-1)
            import matplotlib.pyplot as plt

            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            for i, p in enumerate(pk.clusters):
                print("  Subpeak:", p)
                # if p.monoiso > 6396:
                #     continue
                cplot(p.decon_centroids, color=colors[i % len(colors)], factor=1)
                cplot(p.massdist, color=colors[i % len(colors)], factor=-1)
            plt.show()


    exit()
    # filepath = "C:\\Data\\IsoNN\\test2.txt"
    filepath = "Z:\\Group Share\\JGP\\js8b05641_si_001\\etd_spectrum.txt"
    spectrum = np.loadtxt(filepath, skiprows=0)
    spectrum2 = datachop(spectrum, 891.195, 893.757)


    e1 = eng.encode(spectrum2)
    print("Encoded:", e1.shape)

    # z = eng.predict_charge(spectrum2)
    # print("Predicted Charge:", z)

    # exit()
    pks = eng.process_spectrum(spectrum2)
    print(len(pks.peaks))
    # exit()
    plot_pks(pks, centroids=spectrum, show=True, title="C Interface")
    plt.show()
