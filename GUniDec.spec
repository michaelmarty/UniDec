# -*- mode: python -*-
import os
import pymzml
import datetime
import platform
import zipfile
import time
import fnmatch
from multiprocessing import freeze_support
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules
from os import listdir
from PyInstaller import compat
import matplotlib
import sys
import hashlib
import shutil
import subprocess


PRIVATE_REPOSITORY_ROOT = os.path.abspath(SPECPATH)
PUBLIC_UNIDEC_ROOT = os.path.join(PRIVATE_REPOSITORY_ROOT, 'public', 'UniDec')
PUBLIC_UNIDEC_PACKAGE = os.path.join(PUBLIC_UNIDEC_ROOT, 'unidec')


def public_path(*parts):
    return os.path.join(PUBLIC_UNIDEC_ROOT, *parts)


def hashfile(path):
    # BUF_SIZE is totally arbitrary, change for your app!
    BUF_SIZE = 65536  # lets read stuff in 64kb chunks!

    md5 = hashlib.md5()
    sha256 = hashlib.sha256()

    with open(path, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            md5.update(data)
            sha256.update(data)

    md5hash = md5.hexdigest()
    sha256hash = sha256.hexdigest()

    print("MD5: {0}".format(md5hash))
    print("SHA256: {0}".format(sha256hash))

    return md5hash, sha256hash

freeze_support()


def dir_files(path, rel, type='DATA'):
    ret = []
    for p, d, f in os.walk(path):
        relpath = p.replace(path, '')[1:]
        for fname in f:
            ret.append((os.path.join(rel, relpath, fname),
                        os.path.join(p, fname), type))
    return ret


# Some basics
tstart = time.perf_counter()
system = platform.system()
date = datetime.date.today().strftime("%y%m%d")

# Create names of files and directories
exename = 'GUI_UniDec'
if system == "Windows":
    exename += ".exe"
outputdir = 'UniDec_' + system
zipdirectory = outputdir + "_" + date + ".zip"


hiddenimportslist=['scipy.special._ufuncs_cxx', 'scipy.linalg.cython_blas', 'scipy.linalg.cython_lapack',
                 'scipy.special.cython_special','numpy',
                 'scipy._lib.messagestream','clr', 'clr_loader', 'pythonnet',
                 'encodings', 'encodings.__init__',
                 'packaging', 'packaging.version', 'packaging.specifiers',
                 'pubsub', 'pubsub.core', 'matplotlib.backends.backend_ps', 'matplotlib.backends.backend_pdf',
                 'pycparser', 'numba.core.types.old_scalars']


excludeslist = [
    'IPython', 'statsmodels', 'pyopenms', 'sklearn', 'GdkPixbuf',
    'pyQT4', 'pygobject', 'pygtk', 'pyside', 'PySide2', 'shiboken2',
    'PyQt5', 'torch',
    'pywin32', 'win32api', 'win32con', 'win32gui', 'win32com',
    'win32process', 'win32security', 'win32clipboard',
    'win32service', 'win32event', 'win32file'
]


# Preserve the wheel layouts expected by the native-library loaders:
# isogen/bin/<library>, isodec/bin/<library>, and their binary model files.
native_package_binaries = collect_dynamic_libs('isogen') + collect_dynamic_libs('isodec')
native_package_datas = (
    collect_data_files(
        'isogen',
        includes=['models/*.bin', 'comp_model_parameters/*.txt'],
    )
    + collect_data_files('isodec', includes=['modelparams/*.bin'])
)
importer_hiddenimports = collect_submodules('UniDecImporter')
importer_package_datas = collect_data_files(
    'UniDecImporter',
    includes=[
        'Thermo/*.dll',
        'Thermo/RawFileReaderLicense.doc',
        'Waters/Waters_MassLynxSDK_EULA.txt',
    ],
)


# Analysis of packages
a = Analysis([public_path('unidec', 'Launcher.py')],
             pathex=[PUBLIC_UNIDEC_ROOT],
             binaries=native_package_binaries,
             datas=native_package_datas,
             excludes=excludeslist,
             hiddenimports=hiddenimportslist + importer_hiddenimports,
             hookspath=None,
             runtime_hooks=None)

# UniDecImporter loads its .NET assemblies by absolute path relative to
# UniDecImporter/Thermo/RawFileReader.py. Add them after Analysis so PyInstaller
# preserves that package layout and does not filter them as ordinary binaries.
for source, destination in importer_package_datas:
    bundled_path = os.path.join(destination, os.path.basename(source))
    a.datas.append((bundled_path, source, 'DATA'))



# Add extra things
if system == "Windows":
    a.datas += [('UniDec.exe', public_path('unidec', 'bin', 'UniDec.exe'), 'DATA')]
    a.datas += [('h5repack.exe', public_path('unidec', 'bin', 'h5repack.exe'), 'DATA')]
    a.datas += [('unimod.sqlite', public_path('unidec', 'bin', 'unimod.sqlite'), 'DATA')]
    # a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_def.2.dll', 'BINARY')]
    # a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_avx2.2.dll', 'BINARY')]
    # a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_intel_thread.2.dll', 'BINARY')]
    a.datas += [('pymzml\\version.txt', compat.base_prefix + '\\Lib\\site-packages\\pymzml\\version.txt', 'DATA')]
     # Copy over all the DLLs from the bin folder
    for file in os.listdir(public_path('unidec', 'bin')):
        if fnmatch.fnmatch(file, '*.dll'):
            add = [(file, public_path('unidec', 'bin', file), 'BINARY')]
            a.datas += add
            # print add

    # Waters DLLS
    a.datas += [('Waters_MassLynxSDK_EULA.txt', public_path('unidec', 'bin', 'Waters_MassLynxSDK_EULA.txt'), 'DATA')]
    a.datas += [('CDCReader.exe', public_path('unidec', 'bin', 'CDCReader.exe'), 'DATA')]

#LINUX DEPENDENCIES:
elif system == "Linux":
    a.datas += [('unideclinux', public_path('unidec', 'bin', 'unideclinux'), 'BINARY')]


a.datas += [('cacert.pem', public_path('unidec', 'bin', 'cacert.pem'), 'DATA')]
a.datas += [('logo.ico', public_path('unidec', 'bin', 'logo.ico'), 'DATA')]
a.datas += [('mass_table.csv', public_path('unidec', 'bin', 'mass_table.csv'), 'DATA')]
a.datas += [('metaunidec/images/allButton.png', public_path('unidec', 'metaunidec', 'images', 'allButton.png'), 'DATA')]
a.datas += [('metaunidec/images/peakRightClick.png', public_path('unidec', 'metaunidec', 'images', 'peakRightClick.png'), 'DATA')]
a.datas += [('metaunidec/images/rightClick.png', public_path('unidec', 'metaunidec', 'images', 'rightClick.png'), 'DATA')]
a.datas += [('UniDecLogoMR.png', public_path('unidec', 'bin', 'UniDecLogoMR.png'), 'DATA')]

a.datas.extend(dir_files(os.path.join(os.path.dirname(pymzml.__file__), 'obo'), 'obo'))

a.datas.extend(dir_files(public_path('unidec', 'bin', 'Presets'), 'Presets'))
a.datas.extend(dir_files(public_path('unidec', 'bin', 'Example Data'), 'Example Data'))
a.datas.extend(dir_files(compat.base_prefix + '\\Lib\\site-packages\\matchms\\data', "matchms\\data"))

mkldir = compat.base_prefix + "/Lib/site-packages/numpy/DLLs"
newdir = "numpy/DLLs"
a.datas.extend(dir_files(mkldir, newdir, type="BINARY"))
#a.datas.extend(
#    [(mkldir + "/" + mkl, newdir, 'BINARY') for mkl in listdir(mkldir) if mkl.startswith('mkl_') or mkl.startswith('libio')])
# rdkitlibs = compat.base_prefix + "/Lib/site-packages/rdkit.libs"
# a.datas.extend(dir_files(rdkitlibs, ''))

# Assemble and build
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name=exename,
          debug=True,
          strip=None,
          upx=False,
          console=True, icon=public_path('unidec', 'bin', 'logo.ico'))
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=False,
               name=outputdir)

path = os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist', 'UniDec_Windows', 'GUI_UniDec.exe')
import subprocess


dst = os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist', 'UniDec_Windows', 'obo')
src = os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist', 'UniDec_Windows', '_internal', 'obo')


shutil.copytree(src, dst)
shutil.copy(public_path('readme.md'), os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist', 'UniDec_Windows'))
shutil.copy(public_path('installer.bat'), os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist', 'UniDec_Windows'))


out = subprocess.call(path)
if out != 0:
    exit()

# Pause and ask for yes or no on zipping
userin = input("Is Sign Complete? Zip build? (y/n): ")
if userin.lower() != 'y':
    print("Exiting without zipping.")
    sys.exit(0)

print("Zipping...")
# Zip up the final file
os.chdir(os.path.join(PRIVATE_REPOSITORY_ROOT, 'dist'))

zipf = zipfile.ZipFile(zipdirectory, 'w')
for root, dirs, files in os.walk(outputdir):
    for file in files:
        zipf.write(os.path.join(root, file), compress_type=zipfile.ZIP_DEFLATED)
zipf.close()
print("Zipped to", zipdirectory, "from", outputdir)

hashfile(zipdirectory)

tend = time.perf_counter()
print("Build Time: %.2gm" % ((tend - tstart) / 60.0))
