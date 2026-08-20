# -*- mode: python -*-
import os
import pymzml
import datetime
import platform
import zipfile
import time
import fnmatch
from multiprocessing import freeze_support
from PyInstaller.utils.hooks import collect_submodules
from os import listdir
from PyInstaller import compat
import matplotlib
import sys
import hashlib
import shutil
import subprocess





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


# Analysis of packages
a = Analysis(['unidec\\Launcher.py'],
             pathex=[os.getcwd()],
             excludes=excludeslist,
             hiddenimports=hiddenimportslist,
             hookspath=None,
             runtime_hooks=None)



# Add extra things
if system == "Windows":
    a.datas += [('UniDec.exe', 'unidec\\bin\\UniDec.exe', 'DATA')]
    a.datas += [('CDCReader.exe', 'unidec\\bin\\CDCReader.exe', 'DATA')]
    a.datas += [('h5repack.exe', 'unidec\\bin\\h5repack.exe', 'DATA')]
    a.datas += [('unimod.sqlite', 'unidec\\bin\\unimod.sqlite', 'DATA')]
    #a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_def.2.dll', 'BINARY')]
    #a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_avx2.2.dll', 'BINARY')]
    #a.datas += [('numpy/DLLs', 'unidec\\bin\\mkl_intel_thread.2.dll', 'BINARY')]
    a.datas += [('pymzml\\version.txt', compat.base_prefix + '\\Lib\\site-packages\\pymzml\\version.txt', 'DATA')]
    a.datas += [('massql\\msql.ebnf', compat.base_prefix + '\\Lib\\site-packages\\massql\\msql.ebnf', 'DATA')]

     # Copy over all the DLLs from the bin folder
    for file in os.listdir('unidec\\bin'):
        if fnmatch.fnmatch(file, '*.dll'):
            add = [(file, 'unidec\\bin\\' + file, 'BINARY')]
            a.datas += add
            # print add

    a.datas += [('RawFileReaderLicense.doc', 'unidec\\UniDecImporter\\Thermo\\RawFileReaderLicense.doc', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.Data.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.Data.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.RawFileReader.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.RawFileReader.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.MassPrecisionEstimator.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.MassPrecisionEstimator.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.BackgroundSubtraction.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.BackgroundSubtraction.dll', 'BINARY')]
    a.datas += [('Waters_MassLynxSDK_EULA.txt', 'unidec\\bin\\Waters_MassLynxSDK_EULA.txt', 'DATA')]

    a.datas += [('BaseCommon.dll', 'unidec\\UniDecImporter\\Agilent\\BaseCommon.dll', 'BINARY')]
    a.datas += [('BaseDataAccess.dll', 'unidec\\UniDecImporter\\Agilent\\BaseDataAccess.dll', 'BINARY')]
    a.datas += [('MassSpecDataReader.dll', 'unidec\\UniDecImporter\\Agilent\\MassSpecDataReader.dll', 'BINARY')]
    a.datas += [('BaseTof.dll', 'unidec\\UniDecImporter\\Agilent\\BaseTof.dll', 'BINARY')]
    a.datas += [('BaseError.dll', 'unidec\\UniDecImporter\\Agilent\\BaseTof.dll', 'BINARY')]
    a.datas += [('agtsampleinforw.dll', 'unidec\\UniDecImporter\\Agilent\\agtsampleinforw.dll', 'BINARY')]
    a.datas += [('MIDAC.dll', 'unidec\\UniDecImporter\\Agilent\\MIDAC.dll', 'BINARY')]
    a.datas += [('msvcp120.dll', 'unidec\\UniDecImporter\\Agilent\\msvcp120.dll', 'BINARY')]
    a.datas += [('msvcr120.dll', 'unidec\\UniDecImporter\\Agilent\\msvcr120.dll', 'BINARY')]
    a.datas += [('Interop.shell32.dll', 'unidec\\UniDecImporter\\Agilent\\Interop.shell32.dll', 'BINARY')]
    a.datas += [('BaseDataAccess.dll.config', 'unidec\\UniDecImporter\\Agilent\\BaseDataAccess.dll.config', 'BINARY')]
    a.datas += [('AgtFileSelectionDialog.dll', 'unidec\\UniDecImporter\\Agilent\\AgtFileSelectionDialog.dll', 'BINARY')]

    a.datas += [('BlaisWiff.dll', 'unidec\\UniDecImporter\\Sciex\\BlaisWiff.dll', 'BINARY')]
    a.datas += [('Clearcore2.Compression.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Compression.dll', 'BINARY')]
    a.datas += [('Clearcore2.Data.AnalystDataProvider.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Data.AnalystDataProvider.dll', 'BINARY')]
    a.datas += [('Clearcore2.Data.CommonInterfaces.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Data.CommonInterfaces.dll', 'BINARY')]
    a.datas += [('Clearcore2.Data.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Data.dll', 'BINARY')]
    a.datas += [('Clearcore2.Data.WiffReader.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Data.WiffReader.dll', 'BINARY')]
    a.datas += [('Clearcore2.Muni.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Muni.dll', 'BINARY')]
    a.datas += [('Clearcore2.InternalRawXYProcessing.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.InternalRawXYProcessing.dll', 'BINARY')]
    a.datas += [('Clearcore2.RawXYProcessing.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.RawXYProcessing.dll', 'BINARY')]
    a.datas += [('Clearcore2.StructuredStorage.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.StructuredStorage.dll', 'BINARY')]
    a.datas += [('Clearcore2.Utility.dll', 'unidec\\UniDecImporter\\Sciex\\Clearcore2.Utility.dll', 'BINARY')]
    a.datas += [('Sciex.ClearCore.FMan.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.ClearCore.FMan.dll', 'BINARY')]
    a.datas += [('Sciex.Data.Processing.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.Data.Processing.dll', 'BINARY')]
    a.datas += [('Sciex.Data.SimpleTypes.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.Data.SimpleTypes.dll', 'BINARY')]
    a.datas += [('Sciex.Data.XYData.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.Data.XYData.dll', 'BINARY')]
    a.datas += [('Sciex.FMan.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.FMan.dll', 'BINARY')]
    a.datas += [('Sciex.FMan.UI.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.FMan.UI.dll', 'BINARY')]
    a.datas += [('Sciex.TofTof.T2DFman.dll', 'unidec\\UniDecImporter\\Sciex\\Sciex.TofTof.T2DFman.dll', 'BINARY')]
    a.datas += [('WiffReaderCOM.dll', 'unidec\\UniDecImporter\\Sciex\\WiffReaderCOM.dll', 'BINARY')]
    a.datas += [('zlib.net.dll', 'unidec\\UniDecImporter\\Sciex\\zlib.net.dll', 'BINARY')]



    a.datas += [('RawFileReaderLicense.doc', 'unidec\\UniDecImporter\\Thermo\\RawFileReaderLicense.doc', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.Data.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.Data.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.RawFileReader.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.RawFileReader.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.MassPrecisionEstimator.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.MassPrecisionEstimator.dll', 'BINARY')]
    a.datas += [('ThermoFisher.CommonCore.BackgroundSubtraction.dll', 'unidec\\UniDecImporter\\Thermo\\ThermoFisher.CommonCore.BackgroundSubtraction.dll', 'BINARY')]
    a.datas += [('Waters_MassLynxSDK_EULA.txt', 'unidec\\bin\\Waters_MassLynxSDK_EULA.txt', 'DATA')]


    a.datas += [('isodeclib.dll', 'unidec\\IsoDec\\isodeclib.dll', 'BINARY')]
    a.datas += [('isodeclib.lib', 'unidec\\IsoDec\\isodeclib.lib', 'BINARY')]
    a.datas += [('isogen.dll', 'unidec\\IsoDec\\isogen.dll', 'BINARY')]


    a.datas += [('phase_model_8.bin', 'unidec\\IsoDec\\phase_model_8.bin', 'BINARY')]
    a.datas += [('phase_model_4.bin', 'unidec\\IsoDec\\phase_model_4.bin', 'BINARY')]


    #LINUX DEPENDENCIES:
    a.datas += [('libisodeclib.so', 'C:\\Python\\UniDecDev\\unidec\\IsoDec\\libisodeclib.so', 'BINARY')]
    a.datas += [('libfftw3.so', 'C:\\Python\\UniDecDev\\unidec\\IsoDec\\libfftw3.so', 'BINARY')]
    a.datas += [('libisogen.so', 'C:\\Python\\UniDecDev\\unidec\\IsoDec\\libisogen.so', 'BINARY')]










elif system == "Linux":
    a.datas += [('unideclinux', 'unidec/bin/unideclinux', 'BINARY')]


a.datas += [('cacert.pem', os.path.join('unidec\\bin', 'cacert.pem'), 'DATA')]
a.datas += [('logo.ico', 'unidec\\bin\\logo.ico', 'DATA')]
a.datas += [('mass_table.csv', 'unidec\\bin\\mass_table.csv', 'DATA')]
a.datas += [('metaunidec/images/allButton.png', 'metaunidec\\images\\allButton.png', 'DATA')]
a.datas += [('metaunidec/images/peakRightClick.png', 'metaunidec\\images\\peakRightClick.png', 'DATA')]
a.datas += [('metaunidec/images/rightClick.png', 'metaunidec\\images\\rightClick.png', 'DATA')]
a.datas += [('UniDecLogoMR.png', 'unidec\\bin\\UniDecLogoMR.png', 'DATA')]

a.datas.extend(dir_files(os.path.join(os.path.dirname(pymzml.__file__), 'obo'), 'obo'))

a.datas.extend(dir_files("unidec\\bin\\Presets", 'Presets'))
a.datas.extend(dir_files("unidec\\bin\\Example Data", 'Example Data'))
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
          console=True, icon='unidec\\bin\\logo.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=False,
               name=outputdir)

path = "C:\\Python\\UniDec\\dist\\UniDec_Windows\\GUI_UniDec.exe"
import subprocess


dst = "C:\\Python\\UniDec\\dist\\UniDec_Windows\\obo"
src = "C:\\Python\\UniDec\\dist\\UniDec_Windows\\_internal\\obo"


shutil.copytree(src, dst)
shutil.copy("C:\\Python\\UniDec\\readme.md", "C:\\Python\\UniDec\\dist\\UniDec_Windows")
shutil.copy("C:\\Python\\UniDec\\INSTALLER.bat", "C:\\Python\\UniDec\\dist\\UniDec_Windows")


out = subprocess.call(path)
if out != 0:
    exit()

# exit()

print("Zipping...")
# Zip up the final file
os.chdir("dist")

zipf = zipfile.ZipFile(zipdirectory, 'w')
for root, dirs, files in os.walk(outputdir):
    for file in files:
        zipf.write(os.path.join(root, file), compress_type=zipfile.ZIP_DEFLATED)
zipf.close()
print("Zipped to", zipdirectory, "from", outputdir)

hashfile(zipdirectory)

tend = time.perf_counter()
print("Build Time: %.2gm" % ((tend - tstart) / 60.0))
