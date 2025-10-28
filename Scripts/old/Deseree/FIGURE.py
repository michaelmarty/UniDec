import sys
sys.path.append("C:\\Python\\UniDec")

dir="Z:\Group Share\Deseree\Single belts\\20170524_DJR_POPC_T0mBIOSPIN1_Ramp_0_200_20_1_170524121702\\20170524_DJR_POPC_T0mBIOSPIN1_Ramp_0_200_20_1_170524121702_120_unidecfiles"
file="20170524_DJR_POPC_T0mBIOSPIN1_Ramp_0_200_20_1_170524121702_120_rawdata.txt"

path=os.path.join(dir,file)

dir="Z:\Group Share\Bremen\POPC_D1T0-1m_HCD\\20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file="20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path2=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\Single belts\\20170525_DJR_POPC_T1BIOSPIN1_Ramp_0_200_20_2\\20170525_DJR_POPC_T1BIOSPIN1_Ramp_0_200_20_2_120_unidecfiles"
file="20170525_DJR_POPC_T1BIOSPIN1_Ramp_0_200_20_2_120_rawdata.txt"

path3=os.path.join(dir,file)

dir="Z:\Group Share\Bremen\POPC_D1T1-2m_HCD\\20170207_PB_POPC_ND_D1T1-2m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file="20170207_PB_POPC_ND_D1T1-2m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path4=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\Single Belts for Figure\\201705113_DJR_POPC_T2_Ramp_0_200_20_1\\201705113_DJR_POPC_T2_Ramp_0_200_20_1_140_unidecfiles"
file="201705113_DJR_POPC_T2_Ramp_0_200_20_1_140_rawdata.txt"

path5=os.path.join(dir,file)

dir="Z:\Group Share\Bremen\POPC_D1T0-2m_HCD\\20170207_P1D_POPC_ND_D1T0-2m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file="20170207_P1D_POPC_ND_D1T0-2m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path6=os.path.join(dir,file)


paths=[path, path2, path3,path4, path5,path6]
colors=["black","purple","orange","blue","green","red",]
for i,path in enumerate(paths):
    print(path)

    data=np.loadtxt(path)
    data=datachop(data,11380,11460)
    plt.plot(data[:,0] /1000.0, data[:,1]/np.amax(data[:,1])-i*0.60,color=colors[i])
plt.show()

plt.savefig('Z:\\Group Share\\Deseree\\BeltsShift.pdf')
