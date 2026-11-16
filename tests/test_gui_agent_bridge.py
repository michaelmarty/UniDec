import json
import os
import subprocess
import sys
import tempfile
import textwrap
import threading
import unittest
from pathlib import Path

from tests._test_support import has_gui_display
from unidec.gui_agent_bridge import GuiAgentBridgeClient, GuiAgentBridgeServer


UNIDEC_ROOT = Path(__file__).resolve().parents[1]


class DeferredPresenter:
    def __init__(self):
        self.dispatched = threading.Event()
        self.complete = None

    def dispatch_agent_request(self, action, request, on_complete=None):
        self.complete = on_complete
        self.dispatched.set()


class TestGuiAgentBridgeProtocol(unittest.TestCase):
    def test_overlapping_command_returns_busy(self):
        with tempfile.TemporaryDirectory(prefix="unidec-gui-bridge-") as directory:
            presenter = DeferredPresenter()
            server = GuiAgentBridgeServer(
                presenter,
                Path(directory, "bridge.json"),
                command_timeout=5,
            )
            server.start()
            try:
                client = GuiAgentBridgeClient(server.descriptor_path, timeout=5)
                first_response = {}

                def run_first():
                    first_response.update(client.get_parameters(["startz"]))

                first_thread = threading.Thread(target=run_first)
                first_thread.start()
                self.assertTrue(presenter.dispatched.wait(timeout=2))

                second_response = client.set_parameters({"startz": 5})
                self.assertFalse(second_response["ok"])
                self.assertEqual(second_response["error"]["code"], "busy")

                presenter.complete({"ok": True, "values": {"startz": 1}})
                first_thread.join(timeout=2)
                self.assertFalse(first_thread.is_alive())
                self.assertTrue(first_response["ok"])
            finally:
                server.stop()


@unittest.skipUnless(has_gui_display(), "wxPython requires a graphical display")
class TestLiveGuiAgentBridge(unittest.TestCase):
    def test_visible_parameter_updates_are_atomic_and_full_run_updates_gui(self):
        script = r'''
            import json
            import os
            import shutil
            import tempfile
            import threading
            import time
            from importlib.resources import as_file, files
            from pathlib import Path
            from unittest.mock import patch

            from unidec.GUniDec import UniDecApp
            from unidec.gui_agent_bridge import GuiAgentBridgeClient
            from unidec.modules.unidec_presbase import UniDecPres


            def invoke(app, method, *args):
                outcome = {}

                def call():
                    outcome["response"] = method(*args)

                worker = threading.Thread(target=call)
                worker.start()
                deadline = time.monotonic() + 120
                while worker.is_alive() and time.monotonic() < deadline:
                    app.wx_app.Yield()
                    time.sleep(0.01)
                worker.join(timeout=1)
                assert not worker.is_alive(), "GUI bridge request did not complete"
                return outcome["response"]


            with tempfile.TemporaryDirectory(prefix="unidec-live-gui-bridge-") as directory:
                root = Path(directory)
                descriptor = root / "bridge.json"
                resource = files("unidec").joinpath("bin", "Example Data", "ADH.txt")
                with as_file(resource) as source:
                    spectrum = root / source.name
                    shutil.copy2(source, spectrum)

                with patch.object(UniDecPres, "on_load_default", lambda self, *args, **kwargs: None):
                    app = UniDecApp(ignore_args=True, agent_bridge_file=descriptor)
                try:
                    app.import_config()
                    app.on_open_file(spectrum.name, str(spectrum.parent), clean=True)
                    client = GuiAgentBridgeClient(descriptor, timeout=120)

                    app.view.controls.ctlstartz.SetValue("7")
                    read_response = invoke(app, client.get_parameters, ["startz", "endz"])
                    assert read_response["ok"], read_response
                    assert read_response["revision"] == 0
                    assert read_response["values"]["startz"] == 7
                    assert "startz" in read_response["schema"]

                    update_response = invoke(
                        app,
                        client.set_parameters,
                        {"startz": 10, "endz": 30},
                    )
                    assert update_response["ok"], update_response
                    assert update_response["revision"] == 1
                    assert app.eng.config.startz == 10
                    assert app.eng.config.endz == 30
                    assert app.view.controls.ctlstartz.GetValue() == "10"
                    assert app.view.controls.ctlendz.GetValue() == "30"

                    invalid_response = invoke(
                        app,
                        client.set_parameters,
                        {"startz": 31},
                    )
                    assert not invalid_response["ok"]
                    assert invalid_response["error"]["code"] == "invalid_parameters"
                    assert app.agent_revision == 1
                    assert app.eng.config.startz == 10
                    assert app.eng.config.endz == 30
                    assert app.view.controls.ctlstartz.GetValue() == "10"
                    assert app.view.controls.ctlendz.GetValue() == "30"

                    run_response = invoke(app, client.run, "all")
                    assert run_response["ok"], run_response
                    assert run_response["revision"] == 1
                    assert len(app.eng.data.massdat) > 0
                    assert len(app.eng.pks.peaks) > 0
                    assert app.view.plot1.flag
                    assert app.view.plot2.flag
                    assert app.view.plot4.flag
                    assert app.view.plot6.flag
                    assert app.view.peakpanel.list_ctrl.GetItemCount() == len(app.eng.pks.peaks)
                finally:
                    app.view.Destroy()
                    app.wx_app.Yield()
                    app.wx_app.Destroy()
        '''
        result = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(script)],
            cwd=UNIDEC_ROOT,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
