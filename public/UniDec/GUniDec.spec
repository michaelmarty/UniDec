# -*- mode: python -*-
import datetime
import hashlib
from pathlib import Path
import platform
import shutil
import sys

from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs


root = Path(SPECPATH)
# Match IsoDecGUI's preference for the development checkout.
isodec_checkout = root.parent / "IsoDec"
package_paths = [str(root)]
if (isodec_checkout / "isodec" / "__init__.py").is_file():
    package_paths.insert(0, str(isodec_checkout))
sys.path[:0] = package_paths

system = platform.system()
outputdir = "UniDec_" + system

# Preserve package-relative paths used by the native loaders and importers.
datas = []
binaries = []
for package in ("UniDecImporter", "isodec", "isogen", "pymzml"):
    datas += collect_data_files(package, excludes=["src/**", "tests/**", "**/*.lib"])
    binaries += collect_dynamic_libs(package)

bin_dir = root / "unidec" / "bin"
for name in ("cacert.pem", "logo.ico", "mass_table.csv", "UniDecLogoMR.png",
             "unimod.sqlite", "Waters_MassLynxSDK_EULA.txt"):
    datas.append((str(bin_dir / name), "unidec/bin"))
for name in ("Presets", "Example Data"):
    datas.append((str(bin_dir / name), "unidec/bin/" + name))
datas.append((str(root / "unidec" / "metaunidec" / "images"), "unidec/metaunidec/images"))

if system == "Windows":
    for name in ("UniDec.exe", "CDCReader.exe", "h5repack.exe"):
        binaries.append((str(bin_dir / name), "unidec/bin"))
    binaries += [(str(path), "unidec/bin") for path in bin_dir.glob("*.dll")]
else:
    native_name = "unidecmac" if system == "Darwin" else "unideclinux"
    binaries.append((str(bin_dir / native_name), "unidec/bin"))

hiddenimports = [
    "scipy.special._ufuncs_cxx", "scipy.linalg.cython_blas", "scipy.linalg.cython_lapack",
    "scipy.special.cython_special", "pubsub.core", "matplotlib.backends.backend_ps",
    "matplotlib.backends.backend_pdf", "pycparser",
]
if system == "Windows":
    hiddenimports += ["clr", "clr_loader", "pythonnet"]

a = Analysis(
    [str(root / "unidec" / "Launcher.py")],
    pathex=package_paths,
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["IPython", "statsmodels", "pyopenms", "sklearn", "torch",
              "PyQt5", "PySide2", "shiboken2"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True, name="GUI_UniDec",
    debug=False, strip=False, upx=False, console=True,
    icon=str(bin_dir / "logo.ico"),
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name=outputdir)

# Respect --distpath and package without launching an interactive GUI during builds.
destination = Path(DISTPATH) / outputdir
for name in ("readme.md", "LICENSE", "installer.bat"):
    shutil.copy2(root / name, destination / name)
archive = shutil.make_archive(
    str(Path(DISTPATH) / (outputdir + "_" + datetime.date.today().strftime("%y%m%d"))),
    "zip", root_dir=DISTPATH, base_dir=outputdir,
)
with open(archive, "rb") as stream:
    digest = hashlib.sha256()
    for block in iter(lambda: stream.read(65536), b""):
        digest.update(block)
print("ZIP:", archive)
print("SHA256:", digest.hexdigest())
