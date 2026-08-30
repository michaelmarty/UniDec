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
        self._assert_window_launches("unidec.Launcher", "UniDecLauncher")

    def test_unidec_launches(self):
        self._assert_window_launches("unidec.GUniDec", "UniDecApp")

    def test_unidec_im_launches(self):
        self._assert_window_launches("unidec.UniDecIM", "UniDecIMApp")

    def test_metaunidec_launches(self):
        self._assert_window_launches("unidec.MetaUniDec", "UniDecApp", has_suppression_controls=True)

    def test_unichrom_launches(self):
        self._assert_window_launches("unidec.UniChrom", "ChromApp", has_suppression_controls=True)

    def test_ucd_launches_without_full_stack_button(self):
        self._assert_window_launches("unidec.UniDecCD", "UniDecCDApp", False)

    def test_uccd_launches_with_full_stack_button(self):
        self._assert_window_launches("unidec.UniChromCD", "UniChromCDApp", True)

    def _assert_window_launches(self, module_name, class_name, has_full_stack_button=None,
                                has_suppression_controls=False):
        """Construct a window in an isolated process without entering its event loop."""
        script = f"""
import importlib
from unittest.mock import patch

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
