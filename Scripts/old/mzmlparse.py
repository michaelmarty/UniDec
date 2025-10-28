__author__ = 'Michael.Marty'

import sys
sys.path.append("C:\\Python\\UniDec")
from UniDec.unidec_modules import removeduplicates

timestep=0.25
volts=np.arange(0,100,5)
volts=np.arange(0,100,5)
times=np.arange(0.0,len(volts)*timestep,timestep)
print(volts)
print(times)


directory="Z:\Group Share\Deseree\Montfort Files for Timing"
#directory="C:\\MassLynx\\Mike.PRO\\Data\\150507\\mzML"
#file="20150714_ND_Aqpz_11_POPC_E3_Ramp_0_190_10_2min_01.mzML"
file="20170816_JMD_Montfort_Desalt2_Gas5_interval__Ramp_0_100_5.mzML"
outputheader="2018_Montfort_test_0_100_5"
path=os.path.join(directory,file)

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
            outfile=outputheader+str(volts[v])+".txt"
            outpath=os.path.join(directory,outfile)
            np.savetxt(outpath,data)
            print("Saved:",outpath)
else:
    print("File not found:",path)