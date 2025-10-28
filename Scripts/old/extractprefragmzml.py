__author__ = 'Michael.Marty'

import pymzml
import os

#
#  Change These Value
#

dir="C:\\MassLynx\\Mike.PRO\\Data\\150511\\mzML"

os.chdir(dir)
#files=["MTM_150511_PSPC_50_1_10.mzML"]
files=[f for f in os.listdir(dir) if os.path.isfile(f) and f[-5:]==".mzML"]

#
#  Don't Change These
#

print("Files:",files)
#Loop over all files
for file in files:
    #Setting up input files
    os.chdir(dir)
    path=os.path.join(dir,file)
    print("Reading file:",path)

    #Setting up output files
    header=os.path.splitext(file)[0]
    outfile1=header+"_precursors.txt"
    outfile2=header+"_fragments.txt"
    f1=open(outfile1,"w")
    f2=open(outfile2,"w")

    #Main loop
    msrun=pymzml.run.Reader(path)
    for spectrum in msrun:
        #Get scan time and ms level
        try:
            scantime=float(spectrum["scan time"])
            mslevel=int(spectrum["ms level"])
        except:
            mslevel=0
            scantime=-1

        if mslevel>0 and scantime>=0:
            if mslevel==1:
                #For MS level 1, take five highest peaks and write
                for mz,i in spectrum.highestPeaks(5):
                    sout=str(mz)+","+str(scantime)+","+str(i)+"\n"
                    if i>0:
                        f1.write(sout)
            if mslevel==2:
                #For MS level 2, write all
                for mz,i in spectrum.peaks:
                    sout=str(mz)+","+str(scantime)+","+str(i)+"\n"
                    if i>0:
                        f2.write(sout)

    #Clean Up
    f1.close()
    f2.close()
    print("Done\n")


