import os
import subprocess
import sys
import unittest
from pathlib import Path


UNIDEC_ROOT = Path(__file__).resolve().parents[1]


def _has_gui_display():
    if sys.platform.startswith("win") or sys.platform == "darwin":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


@unittest.skipUnless(_has_gui_display(), "wxPython requires a graphical display")
class TestMajorWindowLaunches(unittest.TestCase):
    def test_launcher_launches(self):
        self._assert_window_launches("unidec.Launcher", "UniDecLauncher", launcher_layout=True)

    def test_unidec_launches(self):
        self._assert_window_launches("unidec.GUniDec", "UniDecApp")

    def test_unidec_im_launches(self):
        self._assert_window_launches("unidec.UniDecIM", "UniDecIMApp")

    def test_metaunidec_launches(self):
        self._assert_window_launches("unidec.MetaUniDec", "UniDecApp", has_suppression_controls=True,
                                     has_chrom_width=False)

    def test_unichrom_launches(self):
        self._assert_window_launches("unidec.UniChrom", "ChromApp", has_suppression_controls=True,
                                     has_chrom_width=True, has_chrom_spectrum_menu=True)

    def test_ucd_launches_without_full_stack_button(self):
        self._assert_window_launches("unidec.UniDecCD", "UniDecCDApp", False)

    def test_uccd_launches_with_full_stack_button(self):
        self._assert_window_launches("unidec.UniChromCD", "UniChromCDApp", True)

    def _assert_window_launches(self, module_name, class_name, has_full_stack_button=None,
                                has_suppression_controls=False, launcher_layout=False,
                                has_chrom_width=None, has_chrom_spectrum_menu=False):
        """Construct a window in an isolated process without entering its event loop."""
        script = f"""
import importlib
import numpy as np
from unittest.mock import Mock, patch

module = importlib.import_module({module_name!r})
app_type = getattr(module, {class_name!r})

# Loading a user's saved default may open a large recent data file. Window
# construction is the behavior under test, so keep the smoke test deterministic.
with patch.object(app_type, "on_load_default", lambda self, *args, **kwargs: None):
    app = app_type(ignore_args=True)

try:
    assert app.view is not None
    expected_full_stack_button = {has_full_stack_button!r}
    if expected_full_stack_button is not None:
        assert hasattr(app.view.controls, "rununidecstack") is expected_full_stack_button
    if {has_suppression_controls!r}:
        suppression_controls = (
            "ctlsuppressiontopn",
            "ctlsuppressiontopx",
            "ctlsuppressionsatellite",
            "ctlsuppressionharmonic",
            "ctlsuppressionstartit",
        )
        assert all(hasattr(app.view.controls, name) for name in suppression_controls)
    expected_chrom_width = {has_chrom_width!r}
    if expected_chrom_width is not None:
        assert hasattr(app.view.controls, "ctldtsig") is expected_chrom_width
        if expected_chrom_width:
            assert app.eng.config.dtsig == 0
            assert hasattr(app.view.controls, "ctlUClineardecon")
            assert hasattr(app.view.controls, "ctlUCtype")
            assert app.view.controls.ctlUClineardecon.GetValue()
            assert app.view.controls.ctlUCtype.GetSelection() == 0
            app.view.controls.ctldtsig.SetValue("2.5")
            app.view.controls.ctlUCtype.SetSelection(1)
            app.view.controls.on_uctype()
            app.view.controls.export_gui_to_config()
            assert app.eng.config.dtsig == 2.5
            assert app.eng.config.UCtype == 1
            assert app.eng.config.UClineardecon == 0
            assert not app.view.controls.ctlUClineardecon.IsEnabled()
    if {launcher_layout!r}:
        buttons = [child for panel in app.view.GetChildren() for child in panel.GetChildren()
                   if child.__class__.__name__ == "Button"]
        labels = [button.GetLabel() for button in buttons]
        assert all("UniDec API Shell" not in label for label in labels)
        im_button = next(button for button in buttons if button.GetLabel().startswith("UniDec IM"))
        assert im_button.GetParent().GetSizer().GetItemPosition(im_button) == (5, 1)
    if {has_chrom_spectrum_menu!r}:
        assert not hasattr(app.view, "open_ud_button")
        assert not hasattr(app.view, "run_ud_button")
        assert not hasattr(app.view, "pick_peaks_button_individual")
        assert not hasattr(app.view, "singlepeakpanel")
        assert not hasattr(app.view, "plot2s")
        assert hasattr(app.view.ypanel, "popupID12")
        plot_sizer = app.view.plotpanel.GetSizer()
        assert plot_sizer.GetItemPosition(app.view.plot2) == (1, 0)
        assert plot_sizer.GetItemSpan(app.view.plot2) == (1, 2)
        assert plot_sizer.GetItemPosition(app.view.plotm) == (2, 0)
        app.view.plot2.plotrefreshtop(np.array([1.0, 2.0]), np.array([1.0, 2.0]), config=app.eng.config)
        assert np.allclose(app.view.plot2.subplot1.get_position().bounds, [0.11, 0.11, 0.86, 0.8])
        raw_data = object()
        spectrum = Mock(rawdata=raw_data)
        app.eng.data.spectra = [spectrum]
        app.eng.filename = "sample.raw"
        app.eng.config.udir = "test-output"
        app.eng.unidec_eng.pass_data_in = Mock()
        launched = Mock()
        with patch("unidec.GUniDec.UniDecApp", return_value=launched) as launch:
            app.on_open_ud(0)
        app.eng.unidec_eng.pass_data_in.assert_called_once_with(
            raw_data, dirname="test-output", fname="sample_spectrum_1.txt")
        launch.assert_called_once()
        launched.start.assert_called_once_with()
finally:
    app.view.Destroy()
    app.wx_app.Yield()
"""

        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=UNIDEC_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )

        self.assertEqual(
            result.returncode,
            0,
            f"{module_name}.{class_name} failed to launch.\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
