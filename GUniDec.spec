# -*- mode: python -*-

import os
import pymzml
import datetime
import sys
import platform
import zipfile
import time
import fnmatch
import comtypes
import encodings
import zipimport
from PyInstaller.utils.hooks import collect_data_files
from multiprocessing import freeze_support

freeze_support()


def dir_files(path, rel):
    ret = []
    for p, d, f in os.walk(path):
        relpath = p.replace(path, '')[1:]
        for fname in f:
            ret.append((os.path.join(rel, relpath, fname),
                        os.path.join(p, fname), 'DATA'))
    return ret


# Some basics
tstart = time.perf_counter()
system = platform.system()
date = datetime.date.today().strftime("%y%m%d")

# Determine if this is a distribution run or internal
if "-nowaters" in sys.argv:
    distmode = True
    print("Distribution Release Mode")
else:
    print("Internal Release Mode")
    distmode = False

# Create names of files and directories
exename = 'GUI_UniDec'
if system == "Windows":
    exename += ".exe"
outputdir = 'UniDec_' + system
if distmode:
    outputdir += '_Dist'
zipdirectory = outputdir + "_" + date + ".zip"

# Analysis of packages
a = Analysis(['Launcher.py'],
             pathex=[os.getcwd()],
             excludes=['pandas', 'IPython', 'Cython', 'statsmodels', 'pyopenms', 'sklearn',
                       'GdkPixbuf', 'PIL', 'pyQT4', 'pygobject', 'pygtk', 'pyside', 'PySide2', 'shiboken2', 'PyQt5'],
             hiddenimports=[  # 'plotly','
                 # 'sklearn', 'sklearn.decomposition', 'sklearn.preprocessing', 'sklearn.utils', 'pytest', 'pluggy',
                 # 'sklearn.utils.testing', 'sklearn.utils._cython_blas',
                 'scipy.special._ufuncs_cxx', 'scipy.linalg.cython_blas', 'scipy.linalg.cython_lapack',
                 'scipy._lib.messagestream',
                 'FileDialog', 'Dialog', 'encodings', 'encodings.__init__',
                 'packaging', 'packaging.version', 'packaging.specifiers',
                 'comtypes', "multiplierz", "comtypes.gen", "comtypes.gen._E7C70870_676C_47EB_A791_D5DA6D31B224_0_1_0",
                 'comtypes.gen._18A9D1D8_42BC_4A5D_AD59_590F4049A5B5_0_8_0',
                 "comtypes.gen.UIAutomationClient", "comtypes.gen.RawReader", "multiplierz.mzAPI.management",
                 # 'wx.lib.pubsub','wx.lib.pubsub.core', 'wx.lib.pubsub.core.kwargs','wx.lib.pubsub.core.publisher',
                 'pubsub', 'pubsub.core.publisherbase', 'pubsub.core.kwargs', 'pubsub.core.kwargs.publisher',
                 'pubsub.core.kwargs.listenerimpl', 'pubsub.core.kwargs.publishermixin',
                 'pubsub.core.listenerbase', 'pubsub.core', 'pubsub.core.kwargs.topicargspecimpl',
                 'pubsub.core.kwargs.topicmgrimpl',
                 'Tkinter', 'FixTk', '_tkinter', 'Tkconstants', 'FileDialog', 'Dialog', 'six',
                 'pymzml.run', 'pymzml.plot', 'pymzml.obo'
                 # , 'requests.packages.chardet.sys', 'requests','urllib3.packages.ordered_dict'
             ],
             hookspath=None,
             runtime_hooks=None)

# Add extra things
if system == "Windows":
    a.datas += [('UniDec.exe', 'unidec_bin\\UniDec.exe', 'DATA')]
    a.datas += [('readme.md', 'readme.md', 'DATA')]
    a.datas += [('LICENSE', 'LICENSE', 'DATA')]
    a.datas += [('libiomp5md.dll', 'unidec_bin\\libiomp5md.dll', 'DATA')]
    a.datas += [('mkl_intel_thread.dll', 'unidec_bin\\mkl_intel_thread.dll', 'DATA')]
    a.datas += [('ucrtbase.dll', 'unidec_bin\\ucrtbase.dll', 'DATA')]
    a.datas += [('vcruntime140.dll', 'unidec_bin\\vcruntime140.dll', 'DATA')]
    a.datas += [('libmypfunc.dll', 'unidec_bin\\libmypfunc.dll', 'DATA')]
    # a.datas += [('rawreader.exe', 'unidec_bin\\rawreader.exe', 'DATA')]
    # a.datas += [('rawreadertim.exe', 'unidec_bin\\rawreadertim.exe', 'DATA')]
    a.datas += [('CDCReader.exe', 'unidec_bin\\CDCReader.exe', 'DATA')]
    a.datas += [('h5repack.exe', 'unidec_bin\\h5repack.exe', 'DATA')]
    a.datas += [('pymzml\\version.txt', 'C:\\Python38\\Lib\\site-packages\\pymzml\\version.txt', 'DATA')]

    for file in os.listdir('unidec_bin'):
        if fnmatch.fnmatch(file, 'api*'):
            add = [(file, 'unidec_bin\\' + file, 'DATA')]
            a.datas += add
            # print add

    if not distmode:
        a.datas += [('MassLynxRaw.dll', 'unidec_bin\\MassLynxRaw.dll', 'DATA')]
        a.datas += [('cdt.dll', 'unidec_bin\\cdt.dll', 'DATA')]
        a.datas += [('Waters_MassLynxSDK_EULA.txt', 'unidec_bin\\Waters_MassLynxSDK_EULA.txt', 'DATA')]
elif system == "Linux":
    a.datas += [('unideclinux', 'unidec_bin/unideclinux', 'DATA')]
    # a.datas += [('unideclinuxIM', 'unidec_bin/unideclinuxIM', 'DATA')]
    a.datas += [('libmypfunc.so', 'unidec_bin/libmypfunc.so', 'DATA')]

a.datas += [('cacert.pem', os.path.join('unidec_bin', 'cacert.pem'), 'DATA')]
a.datas += [('Images/logo.ico', 'logo.ico', 'DATA')]
a.datas += [('metaunidec/logo.ico', 'logo.ico', 'DATA')]
a.datas += [('logo.ico', 'logo.ico', 'DATA')]
a.datas += [('mass_table.csv', 'unidec_bin\\mass_table.csv', 'DATA')]
a.datas += [('metaunidec/images/allButton.png', 'metaunidec\\images\\allButton.png', 'DATA')]
a.datas += [('metaunidec/images/peakRightClick.png', 'metaunidec\\images\\peakRightClick.png', 'DATA')]
a.datas += [('metaunidec/images/rightClick.png', 'metaunidec\\images\\rightClick.png', 'DATA')]
a.datas += [('UniDecLogoMR.png', 'UniDecLogoMR.png', 'DATA')]

a.datas.extend(dir_files(os.path.join(os.path.dirname(pymzml.__file__), 'obo'), 'obo'))

a.datas.extend(dir_files("unidec_bin\\multiplierz", 'multiplierz'))

a.datas.extend(dir_files("unidec_bin\\Presets", 'Presets'))
a.datas.extend(dir_files("unidec_bin\\Example Data", 'Example Data'))

# Can't remember why I needed these...
# grammar=os.path.join(os.path.dirname(lib2to3.__file__),'Grammar.txt')
# a.datas.extend([(os.path.join('lib2to3','Grammar.txt'),grammar,'DATA')])
# grammar=os.path.join(os.path.dirname(lib2to3.__file__),'PatternGrammar.txt')
# a.datas.extend([(os.path.join('lib2to3','PatternGrammar.txt'),grammar,'DATA')])

from os import listdir
from PyInstaller import compat

#mkldir = compat.base_prefix + "/Lib/site-packages/numpy/DLLs"
#a.datas.extend(dir_files(mkldir, ''))
#a.datas.extend(
#    [(mkldir + "/" + mkl, '', 'DATA') for mkl in listdir(mkldir) if mkl.startswith('mkl_') or mkl.startswith('libio')])
#for b in binaries:
#    a.datas += b

# Assemble and build
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name=exename,
          debug=False,
          strip=None,
          upx=False,
          console=True, icon='logo.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=False,
               name=outputdir)

path = "C:\\Python\\UniDec3\\dist\\UniDec_Windows\\GUI_UniDec.exe"
import subprocess

print("Testing Software...", path)

out = subprocess.call(path)

if out != 0:
    exit()

print("Zipping...")
# Zip up the final file
os.chdir("dist")

zipf = zipfile.ZipFile(zipdirectory, 'w')
for root, dirs, files in os.walk(outputdir):
    for file in files:
        zipf.write(os.path.join(root, file), compress_type=zipfile.ZIP_DEFLATED)
zipf.close()
print("Zipped to", zipdirectory, "from", outputdir)

# Final Print
if distmode:
    print("Distribution Release Mode")
else:
    print("Internal Release Mode")

tend = time.perf_counter()
print("Build Time: %.2gm" % ((tend - tstart) / 60.0))
