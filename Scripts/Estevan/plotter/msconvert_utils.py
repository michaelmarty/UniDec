# msconvert_utils.py
from __future__ import annotations
import os
import json
import shutil
import subprocess
import sys
from typing import Optional

CONFIG_PATH = os.path.join(os.path.expanduser("~"), ".plotter_config.json")

def load_config() -> dict:
    try:
        with open(CONFIG_PATH, "r", encoding="utf8") as f:
            return json.load(f)
    except Exception:
        return {}

def save_config(cfg: dict) -> None:
    try:
        with open(CONFIG_PATH, "w", encoding="utf8") as f:
            json.dump(cfg, f, indent=2)
    except Exception:
        pass

def verify_msconvert(path: str) -> bool:
    """
    Quick verification that candidate msconvert is runnable.
    Runs `msconvert --version` and checks output.
    """
    if not path or not os.path.isfile(path):
        return False
    try:
        proc = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=6)
        out = (proc.stdout or "") + (proc.stderr or "")
        return proc.returncode == 0 or "ProteoWizard" in out or "msconvert" in out.lower()
    except Exception:
        return False

def find_msconvert_better(env_path: Optional[str] = None) -> Optional[str]:
    """
    Try multiple strategies to find msconvert.
    Returns full path to msconvert.exe or None.
    """
    # 1) explicit env_path or saved config
    if env_path:
        if verify_msconvert(env_path):
            return env_path

    cfg = load_config()
    saved = cfg.get("msconvert_path")
    if saved and verify_msconvert(saved):
        return saved

    # 2) PROTEOWIZARD_MSCONVERT env var
    env = os.environ.get("PROTEOWIZARD_MSCONVERT")
    if env and verify_msconvert(env):
        cfg["msconvert_path"] = env
        save_config(cfg)
        return env

    # 3) which on PATH
    exe = shutil.which("msconvert")
    if exe and verify_msconvert(exe):
        cfg["msconvert_path"] = exe
        save_config(cfg)
        return exe

    # 4) Windows registry (best-effort, non-fatal)
    if sys.platform.startswith("win"):
        try:
            import winreg
            # look for common keys with InstallLocation or similar
            candidate_keys = [
                r"SOFTWARE\ProteoWizard",
                r"SOFTWARE\WOW6432Node\ProteoWizard",
                r"SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall"
            ]
            for basekey in candidate_keys:
                try:
                    with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, basekey) as k:
                        i = 0
                        while True:
                            try:
                                sub = winreg.EnumKey(k, i)
                            except OSError:
                                break
                            try:
                                with winreg.OpenKey(k, sub) as sk:
                                    try:
                                        val = winreg.QueryValueEx(sk, "InstallLocation")[0]
                                        candidate = os.path.join(val, "msconvert.exe")
                                        if verify_msconvert(candidate):
                                            cfg["msconvert_path"] = candidate
                                            save_config(cfg)
                                            return candidate
                                    except Exception:
                                        pass
                            except Exception:
                                pass
                            i += 1
                except FileNotFoundError:
                    continue
        except Exception:
            pass

    # 5) Common program files locations (shallow search)
    roots = []
    if sys.platform.startswith("win"):
        roots += [os.environ.get("PROGRAMFILES", r"C:\Program Files"),
                  os.environ.get("PROGRAMFILES(X86)", r"C:\Program Files (x86)")]
    else:
        roots += ["/usr/local/bin", "/usr/bin", "/opt"]

    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        # try common path
        candidate = os.path.join(root, "ProteoWizard", "msconvert.exe")
        if verify_msconvert(candidate):
            cfg["msconvert_path"] = candidate
            save_config(cfg)
            return candidate
        # shallow one-level search to avoid long recursion
        try:
            for name in os.listdir(root):
                p = os.path.join(root, name)
                if not os.path.isdir(p):
                    continue
                candidate = os.path.join(p, "msconvert.exe")
                if verify_msconvert(candidate):
                    cfg["msconvert_path"] = candidate
                    save_config(cfg)
                    return candidate
        except Exception:
            continue

    # not found
    return None