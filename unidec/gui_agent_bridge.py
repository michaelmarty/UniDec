"""Authenticated loopback bridge between the MCP process and a live UniDec GUI."""

import hmac
import json
from pathlib import Path
import secrets
import socket
import socketserver
import threading


_MAX_REQUEST_BYTES = 16 * 1024


class GuiAgentBridgeServer:
    """Expose a presenter's allowlisted agent dispatch over localhost."""

    def __init__(self, presenter, descriptor_path, command_timeout=300.0):
        self.presenter = presenter
        self.descriptor_path = Path(descriptor_path).expanduser().resolve()
        self.command_timeout = command_timeout
        self.token = secrets.token_urlsafe(32)
        self._command_lock = threading.Lock()
        self._server = None
        self._thread = None

    def start(self):
        bridge = self

        class Handler(socketserver.StreamRequestHandler):
            def handle(self):
                response = bridge._handle_request(self.rfile.readline(_MAX_REQUEST_BYTES + 1))
                self.wfile.write(json.dumps(response).encode("utf-8") + b"\n")

        class LoopbackServer(socketserver.ThreadingTCPServer):
            allow_reuse_address = True
            daemon_threads = True

        self._server = LoopbackServer(("127.0.0.1", 0), Handler)
        self._thread = threading.Thread(
            target=self._server.serve_forever,
            name="UniDecGuiAgentBridge",
            daemon=True,
        )
        self._thread.start()

        self.descriptor_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor = {
            "version": 2,
            "host": "127.0.0.1",
            "port": self._server.server_address[1],
            "token": self.token,
        }
        self.descriptor_path.write_text(json.dumps(descriptor, indent=2), encoding="utf-8")

    def stop(self):
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
            self._server = None
        try:
            descriptor = json.loads(self.descriptor_path.read_text(encoding="utf-8"))
            if hmac.compare_digest(str(descriptor.get("token", "")), self.token):
                self.descriptor_path.unlink()
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            pass

    def _handle_request(self, payload):
        if not payload or len(payload) > _MAX_REQUEST_BYTES:
            return _error("invalid_request", "The GUI bridge request is empty or too large.")
        try:
            request = json.loads(payload.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            return _error("invalid_request", "The GUI bridge request is not valid JSON.")

        supplied_token = request.get("token")
        if not isinstance(supplied_token, str) or not hmac.compare_digest(supplied_token, self.token):
            return _error("unauthorized", "The GUI bridge token is invalid.")
        action = request.get("action")
        if not isinstance(action, str):
            return _error("invalid_request", "A string action is required.")
        if not self._command_lock.acquire(blocking=False):
            return _error("busy", "The live GUI is already executing an agent command.")

        done = threading.Event()
        outcome = {}

        def complete(response):
            outcome.update(response)
            done.set()

        try:
            self.presenter.dispatch_agent_request(action, request, on_complete=complete)
        except Exception as error:
            self._command_lock.release()
            return _error("dispatch_failed", str(error))

        if not done.wait(self.command_timeout):
            threading.Thread(target=self._release_when_done, args=(done,), daemon=True).start()
            return _error("timeout", "The live GUI command did not finish before the bridge timeout.")

        self._command_lock.release()
        return outcome

    def _release_when_done(self, done):
        done.wait()
        self._command_lock.release()


class GuiAgentBridgeClient:
    """Read a GUI bridge descriptor and invoke one live GUI command."""

    def __init__(self, descriptor_path, timeout=305.0):
        self.descriptor_path = Path(descriptor_path).expanduser().resolve()
        self.timeout = timeout

    def run(self, stage):
        return self._request("run", stage=stage)

    def get_parameters(self, names=None):
        return self._request("get_parameters", names=names)

    def set_parameters(self, values):
        return self._request("set_parameters", values=values)

    def _request(self, action, **values):
        try:
            descriptor = self._read_descriptor()
        except RuntimeError as error:
            return _error("gui_unavailable", str(error))
        payload = {"token": descriptor["token"], "action": action}
        payload.update(values)
        request = json.dumps(payload).encode("utf-8") + b"\n"
        try:
            with socket.create_connection(
                (descriptor["host"], descriptor["port"]), timeout=self.timeout
            ) as connection:
                connection.settimeout(self.timeout)
                connection.sendall(request)
                response = connection.makefile("rb").readline(_MAX_REQUEST_BYTES + 1)
        except OSError as error:
            return _error("gui_unavailable", "Could not connect to the live UniDec GUI: " + str(error))

        if not response or len(response) > _MAX_REQUEST_BYTES:
            return _error("invalid_response", "The live UniDec GUI returned no valid response.")
        try:
            return json.loads(response.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            return _error("invalid_response", "The live UniDec GUI returned invalid JSON.")

    def _read_descriptor(self):
        try:
            descriptor = json.loads(self.descriptor_path.read_text(encoding="utf-8"))
            if descriptor.get("version") != 2:
                raise ValueError("unsupported bridge version")
            host = descriptor["host"]
            port = descriptor["port"]
            token = descriptor["token"]
            if host != "127.0.0.1" or not isinstance(port, int) or not isinstance(token, str):
                raise ValueError("invalid bridge endpoint")
            return descriptor
        except (FileNotFoundError, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise RuntimeError(
                "The live GUI bridge descriptor is missing or invalid: " + str(self.descriptor_path)
            ) from error


def _error(code, message):
    return {"ok": False, "error": {"code": code, "message": message}}
