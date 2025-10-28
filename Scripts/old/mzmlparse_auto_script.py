__author__ = 'Michael.Marty'

import sys
sys.path.append("C:\\Python\\UniDec")
from UniDec.unidec_modules import removeduplicates
import os
import fnmatch


def parse(path,times,timestep,volts,outputheader,directory):
    if os.path.isfile(path):
        for v,time in enumerate(times):
            msrun=pymzml.run.Reader(path)
            data=[]
            i=0
            for spectrum in msrun:
                try:
                    scantime=float(spectrum["scan start time"])
                    #print scantime,time,time+timestep
                except:
                    scantime=-1
                if scantime>=time and scantime<time+timestep:
                    impdat=np.transpose([spectrum.mz,spectrum.i])
                    impdat=impdat[impdat[:,0]>10]
                    data.append(impdat)
                    i+=1
            print("Number of scans:",i," Time: ",time," Length: ",len(data))
            data=np.array(data)
            if len(data)>1:
                try:
                    data=merge_spectra(data)
                except:
                    concat=np.concatenate(data)
                    sort=concat[concat[:,0].argsort()]
                    data=removeduplicates(sort)
                outfile=outputheader+"_"+str(int(volts[v]))+".txt"
                if not os.path.isdir(directory):
                    os.mkdir(directory)
                outpath=os.path.join(directory,outfile)
                np.savetxt(outpath,data)
                print("Saved:",outpath)
    else:
        print("File not found:",path)

def get_files(directory, timestep = 1.0):
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*.mzml'):
            print(file)
            path=os.path.join(directory,file)
            name=os.path.splitext(file)[0]
            newdir=os.path.join(directory, name)
            splits=name.split(sep="_")
            for i,s in enumerate(splits):
                if s.lower() == "ramp":
                    start=float(splits[i+1])
                    end=float(splits[i+2])
                    step=float(splits[i+3])
                    volts=np.arange(start, end+step,step)
                    times = np.arange(0.0, len(volts) * timestep, timestep)
                    parse(path,times,timestep,volts,name,newdir)

directory="Z:\Group Share\Montfort New"
timestep = 1.0
get_files(directory)
exit()


