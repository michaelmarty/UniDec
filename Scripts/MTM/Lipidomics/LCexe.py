import sys, os, subprocess
import LipidListCompare

exe = r"C:\Users\mm96978\AppData\Local\Apps\2.0\ND2QYAB6.QT6\VET5K44V.LYE\skyl..tion_9286511f3362df93_0019.0001_31a2f7d2aa6cc5be\Tools\SkyE3AF\LipidCreator.exe"

topdir = r"Y:\Data\Lipidomics\Libraries"
topdir = r"Z:\Group Share\Annika\Stellar"
os.chdir(topdir)

infile = "LipidList_Comparison_Import.csv"
outfile = "LipidCreatorTransitionList.csv"

inpath = os.path.join(topdir, infile)
outpath = os.path.join(topdir, outfile)

helpcall = [exe, "transitionlist", "help"]
tlcall = [exe, "transitionlist", inpath, outpath, "-h", "0"]
proc = subprocess.Popen(tlcall, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
print(out.decode())
