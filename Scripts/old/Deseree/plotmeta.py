import sys
sys.path.append("C:\\Python\\UniDec")

paths=["Y:\Deseree\Ciara\\20170616_DJR_20_AAControl_POPC_T0m_HCD_Ramp_0_200_20_0_unidecfiles\\20170616_DJR_20_AAControl_POPC_T0m_HCD_Ramp_0_200_20_0_extracts.txt",
       "Y:\Deseree\Ciara\\20170616_DJR_20_IM_POPC_T0m_Ramp_0_200_20_unidecfiles\\20170616_DJR_20_IM_POPC_T0m_Ramp_0_200_20_extracts.txt",
       "Y:\Deseree\Ciara\\20170616_DJR_20_TEA_POPC_T0m_HCD_Ramp_0_200_20_unidecfiles\\20170616_DJR_20_TEA_POPC_T0m_HCD_Ramp_0_200_20_extracts.txt",
       "Y:\Deseree\Ciara\\20170616_DJR_30_IM_Ramp_0_200_20_5_unidecfiles\\20170616_DJR_30_IM_Ramp_0_200_20_5_extracts.txt",
       "Y:\Deseree\Ciara\\20170616_DJR_30_TEA_POPC_T0m_HCD_Ramp_0_200_20_unidecfiles\\20170616_DJR_30_TEA_POPC_T0m_HCD_Ramp_0_200_20_extracts.txt",
       "Y:\Deseree\Ciara\\20170616_DJR_40_AAControl_HCD_Ramp_0_200_20_5_unidecfiles\\20170616_DJR_40_AAControl_HCD_Ramp_0_200_20_5_extracts.txt",

       ]
colors=["blue","blue","blue","red","red","purple",]
xvalues=[0,20,40,60,80,100,120,140,160,180,200]
for i,path in enumerate(paths):
    print(path)
    data=np.loadtxt(path)
    plt.plot(xvalues, data[1:],color=colors[i],marker="o")
plt.show()

#plt.savefig('Z:\\Group Share\\Deseree\\BeltsShift.pdf')
