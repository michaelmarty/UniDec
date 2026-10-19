import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from unidec.modules.CDEng import UniDecCD, cconv2D_preB, fft_fun
from unidec.modules import peakstructure


class TestCDScoring(unittest.TestCase):
    def make_engine(self, charges=(1.0, 2.0)):
        eng = UniDecCD()
        eng.mz = np.arange(1.0, 18.0)
        eng.ztab = np.asarray(charges)
        eng.data.ztab = eng.ztab.copy()
        eng.config.adductmass = 0
        eng.config.massbins = 1
        eng.config.mzbins = 1
        eng.config.CDzbins = 0.5
        eng.config.mzsig = 0
        eng.config.csig = 0
        eng.config.psig = 0
        eng.config.zzsig = 0
        eng.config.msig = 0
        eng.config.beta = 0
        eng.config.numit = 1
        eng.config.rawflag = 1
        eng.config.cdmsflag = 1
        eng.config.datanorm = 0
        eng.X, eng.Y = np.meshgrid(eng.mz, eng.ztab)
        eng.mass = eng.X * eng.Y
        eng.harray = np.zeros_like(eng.X)
        for row, charge in enumerate(eng.ztab):
            # Equal triangular shapes in mass space for each charge.
            eng.harray[row] = np.maximum(1 - np.abs(eng.mass[row] - 8) / 2, 0)
        eng.data.massdat = np.column_stack((np.arange(5., 12.), [0, 0, 1, 2, 1, 0, 0]))
        grid = np.array([0, 0, 0.5, 1, 0.5, 0, 0])[:, None] * np.ones((1, len(charges)))
        eng.data.massgrid = grid.ravel()
        eng.pks = peakstructure.Peaks()
        eng.pks.add_peaks([[8, 2]])
        eng.pks.peaks[0].intervalFWHM = [7, 9]
        eng.pks.peaks[0].badFWHM = False
        eng._set_fit(eng.harray, eng.harray)
        return eng

    def test_ideal_peak_has_unit_components_and_uniscore(self):
        eng = self.make_engine()
        eng.dscore()
        peak = eng.pks.peaks[0]
        for name in ("mscore", "uscore", "cs_score", "fscore", "dscore", "lscore"):
            self.assertAlmostEqual(getattr(peak, name), 1, msg=name)
        self.assertAlmostEqual(eng.pks.uniscore, 1)

    def test_uscore_method_dispatch_and_skip(self):
        eng = self.make_engine()
        self.assertEqual(eng.uscore_method, 1)
        for method in (0, 1, 2):
            with self.subTest(method=method):
                eng.dscore(uscore_method=method)
                self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1)
                self.assertAlmostEqual(eng.pks.uniscore, 1)
        eng.uscore_method = 2
        eng.dscore()
        self.assertAlmostEqual(eng.pks.uniscore, 1)
        with self.assertRaises(ValueError):
            eng.dscore(uscore_method=3)
        eng._clear_scores()
        with patch.object(eng, "make_kernel", side_effect=AssertionError("kernel")):
            eng.pks_uscore(uscore_method=0)
        self.assertEqual(eng.pks.peaks[0].uscore, 1)

    def test_area_method_corrects_ideal_response_without_convolutions(self):
        for shape in (0, 1, 2):
            for rawflag in (0, 1):
                for norm in (0, 1):
                    with self.subTest(shape=shape, rawflag=rawflag, norm=norm):
                        eng = self.make_engine()
                        eng.mz = np.arange(1., 102.)
                        eng.ztab = np.arange(1., 42.)
                        eng.config.CDzbins = 1
                        eng.config.mzsig = 5
                        eng.config.csig = 5
                        eng.config.psfun = shape
                        eng.config.psfunz = shape
                        eng.config.rawflag = rawflag
                        eng.config.datanorm = norm
                        eng.make_kernel()
                        source = np.zeros((41, 101))
                        source[20, 50] = 100
                        observed = np.maximum(cconv2D_preB(source, fft_fun(eng.kernel)), 0)
                        output = (np.maximum(cconv2D_preB(source, fft_fun(eng.mkernel)), 0)
                                  if rawflag == 0 else source.copy())
                        output[output < 1e-12] = 0
                        if norm:
                            output /= output.max()
                        eng.harray = output
                        eng.pks.peaks[0].mass = 51 * 21
                        eng._set_fit(observed, observed)
                        with patch("unidec.modules.CDEng.cconv2D_preB",
                                   side_effect=AssertionError("extra convolution")):
                            eng.pks_uscore(uscore_method=2)
                        self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1, places=8)

    def test_area_method_penalizes_wrong_charge_and_handles_zero_signal(self):
        eng = self.make_engine()
        eng._set_fit(eng.harray[::-1], eng.harray)
        eng.pks_uscore(uscore_method=2)
        self.assertLess(eng.pks.peaks[0].uscore, 1)
        eng.harray[:] = 0
        eng._set_fit(eng.harray, eng.harray)
        eng.pks_uscore(uscore_method=2)
        self.assertEqual(eng.pks.peaks[0].uscore, 0)

    def test_area_method_singleton_fractional_bins_and_clipped_windows(self):
        for charges in ((1.,), (1., 1.5, 2.)):
            with self.subTest(charges=charges):
                eng = self.make_engine(charges)
                eng.pks_uscore(uscore_method=2)
                self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1)
        eng = self.make_engine((1., 2., 3.))
        eng.config.mzsig = 2
        eng.config.csig = 2
        eng.config.psfun = 2
        eng.config.psfunz = 2
        eng.make_kernel()
        source = np.zeros_like(eng.harray)
        source[0, 0] = 1
        observed = np.maximum(cconv2D_preB(source, fft_fun(eng.kernel)), 0)
        eng.harray = source
        eng.pks.peaks[0].mass = 1
        eng._set_fit(observed, observed)
        eng.pks_uscore(uscore_method=2, power=0)
        self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1)

    def test_skipping_uscore_keeps_other_components_and_global_fit_penalty(self):
        eng = self.make_engine()
        eng._set_fit(eng.harray[::-1], eng.harray)
        eng.dscore(uscore_method=0)
        peak = eng.pks.peaks[0]
        self.assertEqual(peak.uscore, 1)
        self.assertAlmostEqual(peak.dscore, peak.mscore * peak.cs_score * peak.fscore)
        self.assertAlmostEqual(eng.pks.uniscore, peak.dscore * eng.config.error)

    def test_charge_resolved_error_is_visible_when_mz_projection_is_identical(self):
        eng = self.make_engine()
        observed = eng.harray[::-1].copy()
        np.testing.assert_allclose(observed.sum(axis=0), eng.harray.sum(axis=0))
        eng._set_fit(observed, eng.harray)
        eng.dscore()
        self.assertLess(eng.pks.peaks[0].uscore, 1)
        self.assertLess(eng.pks.uniscore, 1)

    def test_global_fit_amplitude_is_aligned_before_scoring(self):
        eng = self.make_engine()
        eng._set_fit(eng.harray * 2, eng.harray)
        eng.dscore()
        self.assertAlmostEqual(eng.fit_scale, 2)
        np.testing.assert_allclose(eng.fit2d, eng.observed)
        self.assertAlmostEqual(eng.config.error, 1)
        self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1)
        self.assertAlmostEqual(eng.pks.uniscore, 1)

    def test_singleton_and_fractional_charge_bins(self):
        for charges in ((1.0,), (1.0, 1.5, 2.0)):
            with self.subTest(charges=charges):
                eng = self.make_engine(charges)
                eng.dscore()
                self.assertAlmostEqual(eng.pks.uniscore, 1)

    def test_empty_zero_and_constant_evidence(self):
        eng = self.make_engine()
        for value, expected in ((0., 0.), (2., 1.)):
            evidence = np.full_like(eng.harray, value)
            eng._set_fit(evidence, evidence)
            self.assertEqual(eng.config.error, expected)
        eng.data.massgrid[:] = 0
        eng.harray[:] = 0
        eng._set_fit(eng.harray, eng.harray)
        eng.dscore()
        self.assertEqual(eng.pks.peaks[0].dscore, 0)
        self.assertEqual(eng.pks.uniscore, 0)
        eng.pks = peakstructure.Peaks()
        eng.pks.uniscore = 1
        eng.dscore()
        self.assertEqual(eng.pks.uniscore, 0)

    def test_invalid_evidence_and_missing_fit(self):
        eng = self.make_engine()
        for fitted in (np.zeros((1, 1)), np.full_like(eng.harray, np.nan),
                       np.full_like(eng.harray, np.inf)):
            with self.assertRaises(ValueError):
                eng._set_fit(eng.harray, fitted)
        eng._clear_scores()
        self.assertEqual(eng.fit_scale, 0)
        with self.assertRaisesRegex(ValueError, "deconvolution"):
            eng.dscore()
        self.assertEqual(eng.pks.peaks[0].dscore, 0)

    def test_preserved_evidence_and_scale_invariance(self):
        eng = self.make_engine()
        observed = eng.harray.copy()
        eng._set_fit(observed, eng.harray)
        observed[:] = 0
        self.assertGreater(eng.observed.sum(), 0)
        eng.dscore()
        expected = eng.pks.uniscore
        eng.harray *= 100
        eng._set_fit(eng.observed * 100, eng.fit2d * 100)
        eng.dscore()
        self.assertAlmostEqual(eng.pks.uniscore, expected)

    def test_response_overlap_lowers_uniqueness(self):
        eng = self.make_engine((1.,))
        eng.config.mzsig = 2
        eng.harray[:] = 0
        eng.harray[0, 7] = 1  # mass 8, inside the peak window
        eng.harray[0, 12] = 1  # mass 13, outside it
        eng.make_kernel()
        fit = cconv2D_preB(eng.harray, fft_fun(eng.kernel))
        eng._set_fit(fit, fit)
        eng.pks_uscore()
        self.assertGreater(eng.pks.peaks[0].uscore, 0)
        self.assertLess(eng.pks.peaks[0].uscore, 1)

    def test_profile_response_does_not_convolve_mz_twice_or_change_core_kernel(self):
        eng = self.make_engine((1.,))
        eng.config.mzsig = 2
        eng.config.rawflag = 0
        eng.make_kernel()
        full_kernel = eng.kernel.copy()
        # A narrow mass window includes the complete m/z-profile response.
        eng.harray[0] = np.maximum(1 - np.abs(eng.mz - 8) / 2, 0)
        eng._set_fit(eng.harray, eng.harray)
        eng.pks_uscore()
        self.assertAlmostEqual(eng.pks.peaks[0].uscore, 1)
        np.testing.assert_allclose(eng.kernel, full_kernel)
        self.assertEqual(eng.config.mzsig, 2)

    def test_python_core_preserves_fit_on_measured_scale(self):
        for norm in (0, 1):
            with self.subTest(norm=norm):
                eng = self.make_engine()
                eng.config.datanorm = norm
                eng.harray *= 7
                observed = eng.harray.copy()
                eng.make_kernel()
                eng.decon_core()
                np.testing.assert_allclose(eng.observed, observed)
                np.testing.assert_allclose(eng.fit2d, observed, atol=1e-12)
                self.assertAlmostEqual(eng.config.error, 1)

    def test_native_fit_import_preserves_charge_mz_order(self):
        eng = self.make_engine()
        observed = eng.harray.copy()
        with tempfile.TemporaryDirectory() as directory:
            eng.config.infname = str(Path(directory, "input.txt"))
            eng.config.deconfile = str(Path(directory, "decon.txt"))
            eng.config.fitdatfile = str(Path(directory, "fit.bin"))
            np.savetxt(eng.config.deconfile, observed.T.ravel() * 3)
            observed.T.astype(eng.config.dtype).tofile(eng.config.fitdatfile)
            with patch.object(eng, "export_config"), patch("unidec.tools.unidec_call", return_value=0):
                eng.decon_external_call()
        np.testing.assert_allclose(eng.observed, observed)
        np.testing.assert_allclose(eng.fit2d, observed)
        np.testing.assert_allclose(eng.harray, observed)
        np.testing.assert_allclose(eng.data.fitdat, observed.sum(axis=0))
        eng.dscore()
        self.assertAlmostEqual(eng.pks.uniscore, 1)

    def test_packaged_native_and_python_identity_response_agree(self):
        for norm in (0, 1):
            for rawflag in (0, 1):
                with self.subTest(norm=norm, rawflag=rawflag), tempfile.TemporaryDirectory() as directory:
                    eng = self.make_engine()
                    if not Path(eng.config.UniDecPath).is_file():
                        self.skipTest("Packaged UniDec executable is unavailable")
                    eng.config.datanorm = norm
                    eng.config.rawflag = rawflag
                    eng.harray *= 7
                    observed = eng.harray.copy()
                    eng.config.outfname = str(Path(directory, "native"))
                    eng.config.default_file_names()
                    eng.decon_external_call()
                    native_fit = eng.fit2d.copy()
                    eng.dscore()
                    native_score = eng.pks.uniscore
                    eng.harray = observed
                    eng.make_kernel()
                    eng.decon_core()
                    eng.dscore()
                    np.testing.assert_allclose(eng.fit2d, native_fit, atol=1e-5)
                    self.assertAlmostEqual(eng.pks.uniscore, native_score, places=5)

    def test_packaged_native_nonzero_width_conserves_signal(self):
        with tempfile.TemporaryDirectory() as directory:
            eng = self.make_engine()
            executable = Path(eng.config.UniDecPath)
            if not executable.is_file():
                self.skipTest("Packaged UniDec executable is unavailable")
            native_source = Path(__file__).resolve().parents[1] / "unidec" / "src" / "UniDecCD_Main.c"
            if native_source.stat().st_mtime > executable.stat().st_mtime:
                self.skipTest("Packaged UniDec executable predates the native source change")
            eng.config.mzsig = 2
            eng.config.csig = 1
            eng.config.rawflag = 0
            eng.config.datanorm = 0
            observed_sum = np.sum(eng.harray)
            eng.config.outfname = str(Path(directory, "native"))
            eng.config.default_file_names()

            eng.decon_external_call()

            # Conservation applies to native output before the scoring-only
            # least-squares amplitude adjustment. After one iteration the fit
            # shape is imperfect, so its optimal amplitude need not equal one.
            native_fit = np.fromfile(eng.config.fitdatfile, dtype=eng.config.dtype)
            native_output = np.loadtxt(eng.config.deconfile)
            self.assertAlmostEqual(np.sum(native_fit), observed_sum, places=5)
            self.assertAlmostEqual(np.sum(native_output), observed_sum, places=5)
            native_fit = native_fit.reshape(len(eng.mz), len(eng.ztab)).T
            np.testing.assert_allclose(eng.fit2d, native_fit * eng.fit_scale, rtol=1e-6)
            self.assertAlmostEqual(np.sum(eng.harray), observed_sum, places=5)

    def test_full_python_workflow_scores_and_exports_after_transform(self):
        for rawflag in (0, 1):
            with self.subTest(rawflag=rawflag), tempfile.TemporaryDirectory() as directory:
                eng = self.make_engine()
                eng.exemode = False
                eng.config.rawflag = rawflag
                eng.config.masslb = 0
                eng.config.massub = 40
                eng.config.nativezlb = -1000
                eng.config.nativezub = 1000
                eng.config.peakwindow = 2
                eng.config.peakthresh = 0.1
                eng.config.outfname = str(Path(directory, "workflow"))
                eng.config.default_file_names()
                eng.run_deconvolution(process_data=False)
                eng.pick_peaks()
                self.assertGreater(len(eng.pks.peaks), 0)
                self.assertGreater(eng.pks.uniscore, 0)
                self.assertTrue(all(np.isfinite(p.dscore) for p in eng.pks.peaks))
                eng.export_params()
                export = Path(str(eng.config.outfname) + "_peakparam.dat")
                self.assertIn("DScore", export.read_text().splitlines()[0])
                np.testing.assert_allclose(np.atleast_2d(np.loadtxt(export))[:, -1],
                                           [p.dscore for p in eng.pks.peaks], atol=1e-6)

    def test_native_failure_invalidates_previous_scores(self):
        eng = self.make_engine()
        eng.dscore()
        with tempfile.TemporaryDirectory() as directory:
            eng.config.infname = str(Path(directory, "input.txt"))
            with patch.object(eng, "export_config"), patch("unidec.tools.unidec_call", return_value=9):
                with self.assertRaisesRegex(RuntimeError, "exit code 9"):
                    eng.decon_external_call()
        self.assertIsNone(eng.fit2d)
        self.assertEqual(eng.pks.uniscore, 0)


if __name__ == "__main__":
    unittest.main()
