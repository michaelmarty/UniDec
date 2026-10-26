"""Open UniDecImporter fixtures through the specialized GUI presenters."""

import importlib.util
import json
import os
import platform
import subprocess
import sys
import textwrap
import unittest
from pathlib import Path

from _test_support import has_gui_display, importer_fixture


UNIDEC_ROOT = Path(__file__).resolve().parents[1]
WINDOWS_X64 = platform.system() == "Windows" and platform.machine().lower() in {
    "amd64", "x86_64", "x64"
}


@unittest.skipUnless(has_gui_display(), "wxPython requires a graphical display")
class TestImporterGUIWorkflows(unittest.TestCase):
    def _fixtures(self, specs):
        return [
            (name, str(importer_fixture(self, *relative_parts)))
            for name, relative_parts in specs
        ]

    def _run_gui_script(self, script, fixtures, timeout=300):
        env = os.environ.copy()
        env["UNIDEC_GUI_TEST_FIXTURES"] = json.dumps(fixtures)
        result = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(script)],
            cwd=UNIDEC_ROOT,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_unideccd_opens_supported_importer_formats(self):
        specs = [
            ("txt", ("SingleScan", "test_txt.txt")),
            ("dat", ("SingleScan", "test_dat.dat")),
            ("csv", ("CDMS", "test_csv_cdms.csv")),
            ("npz", ("CDMS", "test_npz_cdms.npz")),
            ("dmt", ("CDMS", "test_dmt_cdms.dmt")),
            ("mzml", ("test_mzml.mzML",)),
            ("mzml-gz", ("test_mzmlgz.mzML.gz",)),
            ("mzxml", ("test_mzxml.mzXML",)),
        ]
        if WINDOWS_X64:
            specs.append(("thermo-raw", ("CDMS", "test_raw_cdms.RAW")))

        self._run_gui_script(
            """
            import json
            import os
            import shutil
            import tempfile
            from pathlib import Path
            from unittest.mock import patch

            import numpy as np
            from unidec.UniDecCD import UniDecCDApp

            def stage(source, root, label):
                source = Path(source)
                folder = root / label
                folder.mkdir()
                target = folder / source.name
                if source.is_dir():
                    shutil.copytree(source, target)
                else:
                    try:
                        os.link(source, target)
                    except OSError:
                        shutil.copy2(source, target)
                return target

            fixtures = json.loads(os.environ["UNIDEC_GUI_TEST_FIXTURES"])
            with patch.object(UniDecCDApp, "on_load_default", lambda self, *args, **kwargs: None):
                app = UniDecCDApp(ignore_args=True)
            try:
                with tempfile.TemporaryDirectory(prefix="unideccd-importers-") as directory:
                    root = Path(directory)
                    for label, source in fixtures:
                        path = stage(source, root, label)
                        with patch.object(app, "on_dataprep_button"):
                            app.on_open_file(None, None, path=str(path), refresh=True)
                        assert app.eng.darray.ndim == 2, label
                        assert app.eng.darray.shape[1] == 5, label
                        assert len(app.eng.farray) > 0, label
                        assert np.isfinite(app.eng.farray).all(), label
                        app.eng.importer.close()
            finally:
                app.view.Destroy()
                app.wx_app.Yield()
                app.wx_app.Destroy()
            """,
            self._fixtures(specs),
        )

    @unittest.skipUnless(
        importlib.util.find_spec("isogen.protein_mods") is not None,
        "IsoDecGUI requires pyisogen with isogen.protein_mods",
    )
    def test_isodec_opens_supported_importer_formats(self):
        specs = [
            ("txt", ("SingleScan", "test_txt.txt")),
            ("dat", ("SingleScan", "test_dat.dat")),
            ("csv", ("SingleScan", "test_csv.csv")),
            ("npz", ("CDMS", "test_npz_cdms.npz")),
            ("dmt", ("CDMS", "test_dmt_cdms.dmt")),
            ("mzml", ("test_mzml.mzML",)),
            ("mzml-gz", ("test_mzmlgz.mzML.gz",)),
            ("mzxml", ("test_mzxml.mzXML",)),
        ]
        if WINDOWS_X64:
            specs.extend([
                ("thermo-raw", ("test_thermo.RAW",)),
                ("waters-raw", ("test_waters.raw",)),
                ("agilent-d", ("test_agilent.d",)),
            ])

        self._run_gui_script(
            """
            import json
            import os
            import sys
            import tempfile
            from pathlib import Path

            import numpy as np
            from UniDecImporter.ImporterFactory import VendorReaderUnavailableError

            sys.argv = ["isodec"]
            from unidec.IsoDecGUI import IsoDecPres

            fixtures = json.loads(os.environ["UNIDEC_GUI_TEST_FIXTURES"])
            app = IsoDecPres()
            try:
                with tempfile.TemporaryDirectory(prefix="isodec-importers-") as directory:
                    for label, source in fixtures:
                        app.eng.outfile = str(Path(directory, label, "output"))
                        try:
                            app.on_open_file(source, refresh=True, silent=True)
                        except VendorReaderUnavailableError as error:
                            if label != "agilent-d":
                                raise
                            print(f"Skipping {label}: {error}")
                            continue
                        finally:
                            app.eng.outfile = None
                        data = app.eng.data.rawdata
                        assert data.ndim == 2, label
                        assert data.shape[1] == 2, label
                        assert len(data) > 0, label
                        assert np.isfinite(data).all(), label
                        app.isodeceng.reader.close()
            finally:
                app.view.Destroy()
                app.wx_app.Yield()
                app.wx_app.Destroy()
            """,
            self._fixtures(specs),
        )

    def test_unichrom_opens_supported_importer_formats(self):
        specs = [
            ("mzml", ("test_mzml.mzML",)),
            ("mzml-gz", ("test_mzmlgz.mzML.gz",)),
            ("mzxml", ("test_mzxml.mzXML",)),
        ]
        if WINDOWS_X64:
            specs.extend([
                ("thermo-raw", ("test_thermo.RAW",)),
                ("waters-raw", ("test_waters.raw",)),
                ("agilent-d", ("test_agilent.d",)),
            ])

        self._run_gui_script(
            """
            import json
            import os
            import shutil
            import tempfile
            from pathlib import Path
            from unittest.mock import patch

            import numpy as np
            from UniDecImporter.ImporterFactory import VendorReaderUnavailableError
            from unidec.UniChrom import ChromApp

            def stage(source, root, label):
                source = Path(source)
                folder = root / label
                folder.mkdir()
                target = folder / source.name
                if source.is_dir():
                    shutil.copytree(source, target)
                else:
                    try:
                        os.link(source, target)
                    except OSError:
                        shutil.copy2(source, target)
                return target

            fixtures = json.loads(os.environ["UNIDEC_GUI_TEST_FIXTURES"])
            with patch.object(ChromApp, "on_load_default", lambda self, *args, **kwargs: None):
                app = ChromApp(ignore_args=True)
            try:
                with tempfile.TemporaryDirectory(prefix="unichrom-importers-") as directory:
                    root = Path(directory)
                    for label, source in fixtures:
                        path = stage(source, root, label)
                        try:
                            app.open_file(str(path), batch=True, refresh=True)
                        except VendorReaderUnavailableError as error:
                            if label != "agilent-d":
                                raise
                            print(f"Skipping {label}: {error}")
                            continue
                        assert app.eng.chromdat.chrom_support, label
                        assert app.eng.ticdat.ndim == 2, label
                        assert app.eng.ticdat.shape[1] == 2, label
                        assert len(app.eng.ticdat) > 0, label
                        assert np.isfinite(app.eng.ticdat).all(), label
                        app.eng.chromdat.close()
            finally:
                app.view.Destroy()
                app.wx_app.Yield()
                app.wx_app.Destroy()
            """,
            self._fixtures(specs),
        )


if __name__ == "__main__":
    unittest.main()
