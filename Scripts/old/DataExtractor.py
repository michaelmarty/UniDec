__author__ = 'Michael.Marty'

import os
import subprocess
import numpy as np

rawreaderpath="C:\\\\cprog\\lib_waters\\MassLynxRaw\\Release\\rawreader.exe"
CDCpath="C:\\\\cprog\\lib_waters\\CDCReader\\Release\\CDCReader.exe"

path="C:\\\\MassLynx\\Mike.PRO\\Data\\270415_stds"
os.chdir(path)
subdirs=[f for f in os.listdir(path) if os.path.isdir(f) and f[-4:]==".raw"]
for fname in subdirs:
    header=os.path.splitext(fname)[0]
    try:
        print(header)
        call=[CDCpath,'-r',fname,'-m',header+"_rawdata.txt",'-i',".\\IMMS\\"+header+"_im.txt",'--mass_start',"750","--mass_end","850",'--ms_bin',"0","--ms_smooth_window","0","--ms_number_smooth","0","--im_bin","0","--sparse","1"]
        subprocess.call(call)
        print("Success: ",fname)
    except:
        print("Failed: ",fname)

exit()


date="180215_SERIES1_"
experiments=['NOCA_','CA_']
samples=['HIS_',"CLIP_"]
pols=['POS_','NEG_']
times=np.arange(0,50,6)

for x in range(0,len(experiments)):
    for i in range(0,len(samples)):
        for t in range(0,len(times)):
            for p in range(0,len(pols)):
                exper=experiments[x]
                samp=samples[i]
                time=times[t]
                header=date+exper+samp+pols[p]+str(time)+"H"
                try:
                    fname=header+'.raw'
                    #fname=os.path.join(path,fname)
                    print(os.path.isdir(fname),fname)
                    call=[CDCpath,'-r',fname,'-m',"test.txt",'-i',".\\IMMS\\"+header+"_im.txt",'--mass_start',"750","--mass_end","850",'--ms_bin',"0","--ms_smooth_window","0","--ms_number_smooth","0","--im_bin","0.1","--sparse","1"]
                    subprocess.call(call)
                    print("Success: ",fname)
                except:
                    print("Failed: ",fname)

exit()
header="Aqpz_P2R2_Ramp5"
functions=list(range(0,20))
volts=list(range(50,250,10))
print(functions)
print(volts)
for i in range(0,len(functions)):
    file=os.path.join(path,header+".raw")
    outfolder=os.path.join(path,header+"_1Dextract")
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    call=[rawreaderpath,"-r",file,"-m",os.path.join(outfolder,header+"_"+str(volts[i])+".txt"),"--ms_bin",str(1),"--fn",str(functions[i])]
    subprocess.call(call)
exit()




path="C:\\\\cprog\\Justin"
header="AA070203"
outputs=[["145",145,170],["170",170,200],["200",200,230],["230",230,260]]
outputs=[["140",140,170],["160",170,200],["180",200,230],["200",230,300]]
print(outputs)

save=np.array(outputs)
np.savetxt(os.path.join(path,header+"slice.csv"),save,fmt="%s",delimiter=",")
for i in range(0,len(outputs)):
    subprocess.call([rawreaderpath,"-r",os.path.join(path,header+".raw"),"-m",os.path.join(os.path.join(path,header+"extract2"),header+"_"+str(outputs[i][0])+".txt"),"--scan_start",str(outputs[i][1]),"--scan_end",str(outputs[i][2]),"--ms_bin",str(1)])
    print(i)



header="FK061214_01"
volt=[30,50,70,90,110,120,140,160,170,180,200]
start=[120,255,450,600,720,1120,1380,1540,2400,1700,2000]
end=[254,449,599,719,1119,1199,1539,1699,2748,1999,2399]
outputs=np.transpose(np.array([volt,start,end]))

header="AA070207"
volt=np.arange(80,205,10)
start=[1,40,80,140,160,200,240,285,341,420,470,520,590]
end=[40,80,110,160,180,225,265,307,420,470,520,570,630]

header="AA070203"
outputs=[["145",145,170],["170",170,200],["200",200,230],["230",230,260]]
header="J300403-6"
outputs=[["4",0,30],["80",30,50],["100",51,70],["120",70,90],[140,91,110],[160,111,130],[180,131,175],[200,176,190]]

header="J080503-4"
outputs=[[4,1,30],[60,170,190],[70,210,230],[80,150,170],[90,230,250],[100,130,150],[110,250,270],[120,110,130],[130,270,290],[140,90,110],[150,290,310],[160,70,90],[170,310,330],[180,50,70],[190,330,350],[200,32,50]]


header="J030204-5"

volt=np.arange(55,205,5)
volt[0]=10
volt[1]=50
volt[2]=60
start=np.arange(10,460,15)
end=np.arange(25,470,15)
start[:6]=start[:6]+5
start[6:]=start[6:]+15
end[6:]=end[6:]+15
end[:6]=end[:6]+5
outputs=np.transpose(np.array([volt,start,end]))